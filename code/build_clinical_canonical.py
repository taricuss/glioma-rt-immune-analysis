#!/usr/bin/env python3
"""
Phase 2.0 adapter — Build clinical_canonical.csv and estimate_scores.csv
for direct use by lmm_interaction_pipeline.py.

Inputs (all verified on disk 2026-09-17):
  1. outputs/intermediate/discovery_clinical_metadata_harmonized.csv
     (patient-level, 35 rows × 7 cols: patient_id, treatment_arm, age, etc.)
  2. data/discovery/GSE326221_series_matrix.txt.gz
     (sample-level GSM ordering + "Patient XX, Timepoint Y" titles +
     care_project_id characteristics line for 75 GSMs)

Outputs:
  1. outputs/intermediate/clinical_canonical.csv — per-sample, matches
     LMM schema: sample_id, patient_id, timepoint, treatment_arm,
     who_grade_initial, codeletion_1p19q, mgmt_methylation,
     estimate_purity, age_at_diagnosis.
  2. outputs/intermediate/estimate_scores.csv — separate score vector
     per sample_id, so LMM can be re-run with custom purity.

Covariate anchoring — REAL sources (2026-10-19 rewrite; pre-spec deviation v1.0):
  The earlier version of this adapter filled 1p19q with a seeded stratified
  shuffle, MGMT with random draws, purity with Beta(12,7) draws, and set every
  grade to "II". Coefficients adjusted for random noise are not defensible, so
  those placeholders are gone. All covariates now come from the Johnson 2026
  Nature Supplementary Tables workbook, joined on
  care_id = <patient_id>T<timepoint_num>  (verified 75/75 exact against
  SupplementaryTable2):

  * timepoint: T1 = primary, T2/T3 = recurrent (per plan §4.3, unchanged).
  * who_grade_initial: SupplementaryTable2.grade of each patient's T1 (initial)
    row, mapped G2/G3/G4 -> II/III/IV, so the covariate is "grade at diagnosis"
    and is constant across a patient's samples.
  * codeletion_1p19q: SupplementaryTable2.subtype. No molecular 1p19q call
    exists in any of the 8 sheets, so Oligodendroglioma -> codel is a
    SUBTYPE-DERIVED PROXY (WHO 2021 defines IDH-mutant + 1p19q-codel as the
    oligodendroglioma integrated diagnosis). Declared as a proxy, not a call.
  * estimate_purity: SupplementaryTable3 malignant nuclei fraction
    (rna_malignant_cell_number / rna_cell_number) — the real tumour-cell
    fraction, per sample/timepoint. Replaces the Beta(12,7) draws.
  * age_at_diagnosis: SupplementaryTable1.case_age_diagnosis_years
    (harmonizer value retained only as a fallback).
  * mgmt_methylation: REMOVED. MGMT appears in none of the 8 sheets, so the
    term is dropped from the locked LMM formula and the model is reported as
    "deviation from pre-spec v1.0" rather than adjusted for a placeholder.

Hard pre-spec 4.5(f) treatment-count lock: in the final canonical CSV,
after excluding P11+P85 (NA arm), we MUST end up with:
  Any_RT = 20 unique patients  (exactly)
  No_RT  = 13 unique patients  (exactly)
If not, abort with exit code 5 (per 4.5(f) rule — update tokens, NEVER
change expected counts).
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "outputs" / "intermediate"
OUT.mkdir(parents=True, exist_ok=True)

CSV_IN_PATIENT = OUT / "discovery_clinical_metadata_harmonized.csv"
SERIES_MATRIX = PROJECT / "data" / "discovery" / "GSE326221_series_matrix.txt.gz"
SUPP_XLSX = PROJECT / "data" / "discovery" / "41586_2026_10612_MOESM3_ESM.xlsx"
EXCLUDED_JSON = OUT / "excluded_patients_treatment_na.json"

OUT_CANON = OUT / "clinical_canonical.csv"
OUT_PURITY = OUT / "estimate_scores.csv"
OUT_AUDIT = OUT / "clinical_canonical_adapter_audit.json"

# Plan §4.3 lock
EXPECTED_ANY_RT_N = 20
EXPECTED_NO_RT_N = 13
EXPECTED_EXCLUDED_N = 2
EXPECTED_UNIQUE_PATIENTS_IN_CANON = 33  # = 35 total - 2 excluded NA

# Covariate anchors
N_CODELETED_IN_33_COMPLETE = 9  # per Johnson 2026 + clinical harmonizer
#                                Any_RT_IDH_complete=True (both ≥5)

T1_IS_PRIMARY = True  # plan sec 4.3: T1 = initial diagnosis, T2/T3 = recurrence


def _parse_series_matrix(series_gz: Path) -> pd.DataFrame:
    """Return one row per GSM with: gsm, patient_id (title), timepoint (title)."""
    with gzip.open(series_gz, mode="rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    gsm_list: list[str] | None = None
    title_list: list[str] | None = None
    care_list: list[str] | None = None

    for line in lines:
        if not line.startswith("!Sample_"):
            continue
        stripped = line.rstrip("\n")
        key = stripped.split("\t", 1)[0].strip()
        # Tokenize quoted TSV row (sample metadata columns = quoted, tab-separated)
        tokens = re.findall(r'"([^"]*)"', stripped)
        if key == "!Sample_geo_accession":
            gsm_list = tokens
        elif key == "!Sample_title":
            title_list = tokens
        # "care project_id: P##T#" lives inside characteristics_ch1
        # characteristics_ch1 repeats many times; the last field of each is "care project_id: X"
        elif key == "!Sample_characteristics_ch1":
            if tokens and all("care project_id" in t for t in tokens):
                care_list = [t.split(":", 1)[1].strip() for t in tokens]

    if gsm_list is None or title_list is None:
        raise RuntimeError("Failed to find GSM or Sample_title rows in series_matrix")

    n = len(gsm_list)
    if len(title_list) != n:
        raise RuntimeError(f"GSM list ({n}) != Sample_title list ({len(title_list)})")
    if care_list is not None and len(care_list) != n:
        raise RuntimeError(f"care_project_id list ({len(care_list)}) != GSM list ({n})")

    rows = []
    re_tp = re.compile(r"^Patient\s+([A-Za-z]?\d+),\s*Timepoint\s*(\d+)")
    # MDA/N-prefix patients have only care_project_id titles (not "Patient XX");
    # care list will be authoritative for them instead.
    for i, gsm in enumerate(gsm_list):
        title = title_list[i]
        care = care_list[i] if care_list is not None else ""
        m = re_tp.match(title.strip())
        if m:
            patient_from_title = f"P{m.group(1)}"
            tp_num = int(m.group(2))
        else:
            patient_from_title = ""
            tp_num = None
        # care_project_id: always last 1-2 digit = timepoint
        m2 = re.match(r"^([A-Z]*\d+?)T(\d+)$", care)
        if m2:
            patient_from_care = m2.group(1)
            tp_from_care = int(m2.group(2))
        else:
            patient_from_care = ""
            tp_from_care = None
        # Priority: care_project_id (authoritative because it's in filename
        # regexes too for MDA/N prefixes), fall back to title
        patient_id = patient_from_care or patient_from_title
        timepoint = tp_from_care if tp_from_care is not None else tp_num
        if not patient_id or timepoint is None:
            raise RuntimeError(
                f"Could not determine patient/timepoint for GSM {gsm}: "
                f"title={title!r}, care_project_id={care!r}"
            )
        rows.append({
            "sample_id": gsm,
            "patient_id": patient_id,
            "timepoint_num": int(timepoint),
        })
    df = pd.DataFrame(rows)
    assert df["sample_id"].nunique() == n, "Duplicate GSMs parsed"
    return df


def main() -> int:
    if not CSV_IN_PATIENT.exists():
        print(f"[!!] Missing patient CSV: {CSV_IN_PATIENT}. Run acquire_clinical_metadata.py first.")
        return 2
    if not SERIES_MATRIX.exists():
        print(f"[!!] Missing series_matrix: {SERIES_MATRIX}")
        return 2

    patient_df = pd.read_csv(CSV_IN_PATIENT)
    if EXCLUDED_JSON.exists():
        _excl_raw = json.loads(EXCLUDED_JSON.read_text(encoding="utf-8"))
        if isinstance(_excl_raw, dict) and "patients" in _excl_raw:
            excluded: list[str] = [p["patient_id"] for p in _excl_raw["patients"]]
        elif isinstance(_excl_raw, list):
            excluded = list(_excl_raw)
        else:
            excluded = []
    else:
        excluded = []
        print(f"[W] {EXCLUDED_JSON} not found; assuming empty exclusion list.")

    gsm_df = _parse_series_matrix(SERIES_MATRIX)
    print(f"[OK] Parsed {len(gsm_df)} GSM sample rows from series_matrix.")

    # T1 -> primary, T2/T3/... -> recurrent (plan sec 4.3)
    gsm_df["timepoint"] = np.where(gsm_df["timepoint_num"] == 1, "primary", "recurrent")

    # --------- Treatment join + exclusions (plan §4.3 lock) ----------
    canon = gsm_df.merge(patient_df, on="patient_id", how="left")
    if canon["treatment_arm"].isna().any():
        missing_patients = sorted(canon.loc[canon["treatment_arm"].isna(), "patient_id"].unique().tolist())
        extra = [p for p in missing_patients if p not in excluded]
        if extra:
            raise RuntimeError(
                f"GSM rows with UNEXPECTED missing treatment_arm patient: {extra}. "
                f"(Allowed missing = {excluded}; found extra outside excluded list). "
                f"Check patient_id parsing (MDA/N/P prefixes all OK?)."
            )
        print(
            f"[OK] {missing_patients} have treatment_arm=NaN in harmonizer CSV — exactly the "
            f"§4.5(f) excluded set. Will be dropped by the exclusion filter next step."
        )

    # Exclude NA-treatment-arm patients (plan §4.3: 2 patients P11+P85)
    canon = canon.loc[~canon["patient_id"].isin(excluded)].copy()
    # Also explicit filter on arm != NA just in case
    canon = canon.loc[canon["treatment_arm"] != "NA"].copy()

    # ---------- 4.5(f) treatment count HARD ABORT --------------------
    uniq_pat = canon.drop_duplicates(subset=["patient_id"])
    n_any_rt = int((uniq_pat["treatment_arm"] == "Any_RT").sum())
    n_no_rt  = int((uniq_pat["treatment_arm"] == "No_RT").sum())
    n_excluded = len(excluded)
    n_uniq_canon = int(uniq_pat.shape[0])
    print(f"[OK] Post-exclusion unique patients: {n_uniq_canon}")
    print(f"[OK]   Any_RT = {n_any_rt}, No_RT = {n_no_rt}, excluded NA = {n_excluded}")
    mis = []
    if n_any_rt != EXPECTED_ANY_RT_N:
        mis.append(f"Any_RT: found {n_any_rt}, expected {EXPECTED_ANY_RT_N}")
    if n_no_rt != EXPECTED_NO_RT_N:
        mis.append(f"No_RT: found {n_no_rt}, expected {EXPECTED_NO_RT_N}")
    if n_excluded != EXPECTED_EXCLUDED_N:
        mis.append(f"excluded NA: found {n_excluded}, expected {EXPECTED_EXCLUDED_N}")
    if n_uniq_canon != EXPECTED_UNIQUE_PATIENTS_IN_CANON:
        mis.append(f"total unique patients in CSV: found {n_uniq_canon}, expected {EXPECTED_UNIQUE_PATIENTS_IN_CANON}")
    if mis:
        print("\n[!!] Pre-spec §4.5(f) treatment-count lock FAILURE. Fix patient_id parsing or update")
        print("     treatment classifier tokens. NEVER change expected counts.")
        for m in mis:
            print(f"     - {m}")
        print("[!!] HARD ABORT (exit 5, per 4.5(f) rule).")
        return 5
    print("[OK] §4.5(f) 20/13/2/33 count lock PASSED.")

    # ---------- Covariate anchoring (REAL supp-table sources) ---------
    # 2026-10-19: all seeded placeholder fills removed. Join key = care_id,
    # verified 75/75 exact against SupplementaryTable2.
    if not SUPP_XLSX.exists():
        print(f"[!!] Missing Johnson supp tables: {SUPP_XLSX}")
        return 2
    supp1 = pd.read_excel(SUPP_XLSX, "SupplementaryTable1", header=1)
    supp2 = pd.read_excel(SUPP_XLSX, "SupplementaryTable2", header=1)
    supp3 = pd.read_excel(SUPP_XLSX, "SupplementaryTable3", header=1)
    supp1["patient_id"] = supp1["patient_id"].astype(str).str.strip()
    for _df in (supp2, supp3):
        _df["care_id"] = _df["care_id"].astype(str).str.strip()

    canon["care_id"] = (canon["patient_id"].astype(str) + "T"
                        + canon["timepoint_num"].astype(int).astype(str))

    GRADE_MAP = {"G2": "II", "G3": "III", "G4": "IV"}

    # (1) grade + subtype per sample (SupplementaryTable2)
    canon = canon.merge(
        supp2[["care_id", "subtype", "grade"]].drop_duplicates("care_id"),
        on="care_id", how="left", validate="many_to_one",
    )
    n_missing_grade = int(canon["grade"].isna().sum())
    n_missing_subtype = int(canon["subtype"].isna().sum())

    # who_grade_initial = each patient's T1 (initial-diagnosis) grade, constant
    # across that patient's samples.
    t1_grades = (canon.loc[canon["timepoint_num"] == 1, ["patient_id", "grade"]]
                 .drop_duplicates("patient_id").set_index("patient_id")["grade"])
    canon["who_grade_initial"] = canon["patient_id"].map(t1_grades).map(GRADE_MAP)
    n_missing_initial_grade = int(canon["who_grade_initial"].isna().sum())

    # (2) codeletion_1p19q: SUBTYPE-DERIVED PROXY (no molecular 1p19q call exists).
    canon["codeletion_1p19q"] = np.where(
        canon["subtype"].astype(str).str.strip() == "Oligodendroglioma",
        "codel", "non_codel",
    )

    # (3) estimate_purity: real malignant-nuclei fraction (SupplementaryTable3),
    # per sample/timepoint.
    _s3 = supp3[["care_id", "rna_malignant_cell_number", "rna_cell_number"]].copy()
    _s3["estimate_purity"] = (
        pd.to_numeric(_s3["rna_malignant_cell_number"], errors="coerce")
        / pd.to_numeric(_s3["rna_cell_number"], errors="coerce")
    )
    canon = canon.merge(_s3[["care_id", "estimate_purity"]],
                        on="care_id", how="left", validate="many_to_one")
    n_missing_purity = int(canon["estimate_purity"].isna().sum())

    # (4) age_at_diagnosis: SupplementaryTable1 authoritative, harmonizer fallback.
    age_map = supp1.set_index("patient_id")["case_age_diagnosis_years"]
    canon["age_at_diagnosis"] = (
        pd.to_numeric(canon["patient_id"].map(age_map), errors="coerce")
        .fillna(pd.to_numeric(canon["age_at_diagnosis"], errors="coerce"))
    )
    n_missing_age = int(canon["age_at_diagnosis"].isna().sum())

    # (5) mgmt_methylation: intentionally NOT created — absent from all supp sheets.
    if n_missing_grade or n_missing_subtype or n_missing_initial_grade \
            or n_missing_purity or n_missing_age:
        print(f"[!!] Covariate join produced missing values — grade={n_missing_grade}, "
              f"subtype={n_missing_subtype}, initial_grade={n_missing_initial_grade}, "
              f"purity={n_missing_purity}, age={n_missing_age}")
        return 6
    print("[OK] Real covariates joined 0 missing: grade/subtype/purity/age from SuppTables 1-3.")
    print(f"       who_grade_initial dist: {canon['who_grade_initial'].value_counts().to_dict()}")
    print(f"       codeletion_1p19q dist : {canon['codeletion_1p19q'].value_counts().to_dict()}")
    print(f"       purity median={canon['estimate_purity'].median():.3f} "
          f"min={canon['estimate_purity'].min():.3f}")
    print("       mgmt_methylation: DROPPED (not in any supp table) -> model = pre-spec deviation v1.0")

    # --------- Write output 1: clinical_canonical.csv ----------------
    keep_cols = [
        "sample_id", "patient_id", "timepoint", "treatment_arm",
        "who_grade_initial", "age_at_diagnosis", "codeletion_1p19q",
        "estimate_purity",
    ]
    canon2 = canon.loc[:, keep_cols].copy()
    canon2 = canon2.sort_values(["patient_id", "sample_id"], kind="mergesort").reset_index(drop=True)
    canon2.to_csv(OUT_CANON, index=False)
    print(f"[OK] Wrote {OUT_CANON}  ({canon2.shape[0]} samples × {canon2.shape[1]} cols)")
    print(f"       Timepoints: {canon2['timepoint'].value_counts().to_dict()}")
    print(f"       Arms × timepoint cross-tab:")
    xtab = pd.crosstab(canon2.timepoint, canon2.treatment_arm)
    with pd.option_context("display.width", 120):
        print(xtab.to_string())

    # --------- Write output 2: estimate_scores.csv ------------------
    purity_df = canon2[["sample_id", "estimate_purity"]].rename(
        columns={"estimate_purity": "tumor_purity"}
    )
    # Add canonical ESTIMATE companion columns (immune / stromal not needed by LMM but
    # keeps schema aligned with estimate_purity.py output if we re-run it later)
    purity_df["immune_score"] = np.nan
    purity_df["stromal_score"] = np.nan
    purity_df["estimate_score"] = np.nan
    purity_df.to_csv(OUT_PURITY, index=False)
    print(f"[OK] Wrote {OUT_PURITY} ({len(purity_df)} rows)")
    print(f"       Purity stats: median={purity_df['tumor_purity'].median():.3f}, "
          f"q1={purity_df['tumor_purity'].quantile(0.25):.3f}, "
          f"q3={purity_df['tumor_purity'].quantile(0.75):.3f}")

    # --------- Write audit ------------------------------------------
    grade_dist = canon2["who_grade_initial"].value_counts(dropna=False).to_dict()
    codel_dist = canon2["codeletion_1p19q"].value_counts(dropna=False).to_dict()
    audit = {
        "ran_at": pd.Timestamp.utcnow().isoformat(timespec="seconds") + "Z",
        "sec4p5f_count_lock_passed": True,
        "covariate_source": "Johnson 2026 Nature Supp Tables 1-3 "
                            "(data/discovery/41586_2026_10612_MOESM3_ESM.xlsx)",
        "covariate_join_key": "care_id = <patient_id>T<timepoint_num>; matched 75/75 vs SuppTable2",
        "covariate_placeholders_removed": [
            "codeletion_1p19q seeded stratified shuffle (seed 7331)",
            "mgmt_methylation random draws at assumed rates (seed 20241101)",
            "estimate_purity Beta(12,7) draws (seed 2901)",
            "who_grade_initial hardcoded 'II'",
        ],
        "prespec_deviation": (
            "v1.0-deviation: mgmt_methylation term dropped from the locked LMM "
            "formula because MGMT is absent from all 8 supplementary sheets."
        ),
        "who_grade_initial": {
            "source": "SupplementaryTable2.grade, each patient's T1 (initial) row",
            "map": GRADE_MAP,
            "distribution": grade_dist,
        },
        "codeletion_1p19q": {
            "source": "SupplementaryTable2.subtype used as a PROXY "
                      "(no molecular 1p19q call available in any sheet)",
            "map": {"Oligodendroglioma": "codel", "Astrocytoma": "non_codel"},
            "distribution": codel_dist,
        },
        "mgmt_methylation": "DROPPED - not present in any of the 8 supplementary sheets.",
        "estimate_purity": {
            "source": "SupplementaryTable3 rna_malignant_cell_number / rna_cell_number",
            "median": round(float(canon2["estimate_purity"].median()), 4),
            "q1": round(float(canon2["estimate_purity"].quantile(0.25)), 4),
            "q3": round(float(canon2["estimate_purity"].quantile(0.75)), 4),
            "n_below_0.05": int((canon2["estimate_purity"] < 0.05).sum()),
        },
        "age_at_diagnosis": {
            "source": "SupplementaryTable1.case_age_diagnosis_years (harmonizer fallback)",
        },
        "treatment_arm_source": "discovery_clinical_metadata_harmonized.csv (unchanged)",
        "counts": {
            "n_samples_total": int(canon2.shape[0]),
            "n_patients_unique": int(uniq_pat.shape[0]),
            "n_any_rt_patients": n_any_rt,
            "n_no_rt_patients": n_no_rt,
            "n_excluded_na_treatment": n_excluded,
        },
        "excluded_patients": excluded,
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, default=str)
    print(f"[OK] Wrote audit -> {OUT_AUDIT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
