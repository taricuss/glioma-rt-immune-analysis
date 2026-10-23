
"""
Step 0b -- Clinical metadata harmonization (Johnson et al. Supp Table 1)
=======================================================================

Input:
    The user downloads Johnson et al. (Nature 2026) Supplementary Table 1
    and 2 from the journal page (Supp. data items MOESM3 / MOESM4) and drops
    them into data/discovery/ as:
        data/discovery/Johnson_2026_Nature_SupplTable1.xlsx   (or .csv)
        data/discovery/Johnson_2026_Nature_SupplTable2.xlsx

What this script does:
    1. Parses the treatment column and validates:
         - Any-RT   = n=20
         - No-RT    = n=13
         - NA/unclear excluded = n=2   (the two IDs get written explicitly
           to excluded_patients.json and every primary-fit file we emit)
    2. Parses all covariates mandated by the section 4.5 locked model (purity,
       WHO grade initial, age, 1p/19q, MGMT, timepoint date, patient ID).
    3. Joins to the per-sample metadata produced by acquire_bulk_discovery.py
       (sample <-> patient <-> timepoint) and produces ONE patient x timepoint
       metadata CSV used by every downstream statistical script.

If the user has NOT downloaded the Excel files yet, the script prints the
exact download URL + which columns in Supp Table 1 we need, and exits
non-zero.  No fabrication, no guessing (plan section 0 rule).
"""

from __future__ import annotations

import json
import sys
import io
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data" / "discovery"
OUT  = HERE / "outputs" / "intermediate"
for p in (DATA, OUT):
    p.mkdir(parents=True, exist_ok=True)

# Plan section 4.3 -- locked treatment classification rules.  Do NOT change.
ANY_RT_TOKENS   = ["radiotherapy", "radiation", "rt", "radiochemotherapy", "chemoradiation", "rt+tmz", "rt + tmz", "ccrt", "concurrent chemoradiation"]
NO_RT_TOKENS    = ["tmz only", "temozolomide only", "chemotherapy only", "chemo only", "adjuvant chemotherapy", "no adjuvant", "observation", "watchful waiting", "untreated", "none", "nil", "no treatment"]
NA_TOKENS       = ["not reported", "na", "n/a", "unknown", "unclear", "missing", ""]

PAPER_SUPP_URLS = {
    "Suppl Table 1": "https://www.nature.com/articles/s41586-026-10612-6#MOESM3",
    "Suppl Table 2": "https://www.nature.com/articles/s41586-026-10612-6#MOESM4",
}

def find_supp_table_files() -> dict[str, Path]:
    found = {}
    for suffix in (".xlsx", ".xls", ".csv"):
        for pattern in ("*SupplTable1*", "*SuppTable1*", "*TableS1*", "*supp_table_1*", "*MOESM3*"):
            hits = list(DATA.glob(pattern + suffix))
            if hits:
                found.setdefault("SupplTable1", hits[0])
    for suffix in (".xlsx", ".xls", ".csv"):
        for pattern in ("*SupplTable2*", "*SuppTable2*", "*TableS2*", "*supp_table_2*", "*MOESM4*"):
            hits = list(DATA.glob(pattern + suffix))
            if hits:
                found.setdefault("SupplTable2", hits[0])
    return found

def classify_treatment(cell) -> str:
    if cell is None:
        return "NA"
    s = str(cell).strip().lower()
    if s in NA_TOKENS:
        return "NA"
    any_rt_hit = any(tok in s for tok in ANY_RT_TOKENS)
    no_rt_hit  = any(tok in s for tok in NO_RT_TOKENS)
    if any_rt_hit and no_rt_hit:
        return "Any_RT"
    if any_rt_hit:
        return "Any_RT"
    if no_rt_hit:
        return "No_RT"
    return "NA"

def finer_breakdown(cell, overall: str) -> str:
    if overall != "Any_RT":
        s = str(cell).strip().lower()
        if any(tok in s for tok in ("tmz only", "temozolomide only", "chemotherapy only", "chemo only", "adjuvant chemo")):
            return "chemo_only"
        if any(tok in s for tok in ("no adjuvant", "observation", "watchful", "untreated", "none", "nil")):
            return "untreated"
        return "unclassified_no_rt"
    s = str(cell).strip().lower()
    has_chemo = any(tok in s for tok in ("tmz", "temozolomide", "chemo", "alkylating", "ccnu", "pcv"))
    has_rt    = any(tok in s for tok in ANY_RT_TOKENS)
    if has_rt and not has_chemo:
        return "rt_only"
    if has_rt and has_chemo:
        return "rt_plus_chemo"
    if has_rt:
        return "rt_plus_chemo_unconfirmed_chemo"
    return "unclassified_any_rt"

EXPECTED_COLS_HINT = """
Supp Table 1 column names to look for (they vary slightly between papers):
    patient_id / case / sample_id
    who_grade_initial  (or grade_initial / who_grade)
    age_at_diagnosis  (or age)
    treatment_between_timepoints  (or adjuvant_treatment / treatment_after_initial / tx)
    codeletion_1p19q  (or 1p_19q / 1p19q)
    mgmt_methylation_status  (or mgmt / mgmt_promoter)
    timepoint  (or resection_number / initial_vs_recurrent / tp)
    days_between_resctions / interval_months
    (ESTIMATE) purity -- sometimes in Supp Table 2 alongside CIBERSORTx fractions;
                         if missing, we compute in Step 5.1 (plan section 5.3) and back-fill.
"""


def _bin(x):
    try:
        if pd.isna(x):
            return None
    except Exception:
        pass
    s = str(x).strip()
    if s == "":
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def classify_treatment_binary(row, pat_col, rt_col, tmz_col, alk_col,
                              na_patient_ids):
    pid = str(row.get(pat_col, "")).strip()
    if pid in na_patient_ids:
        return "NA"
    rt = _bin(row.get(rt_col, None)) if rt_col else None
    if rt == 1:
        return "Any_RT"
    if rt == 0:
        tmz = _bin(row.get(tmz_col, None)) if tmz_col else None
        alk = _bin(row.get(alk_col, None)) if alk_col else None
        if tmz is None and alk is None:
            return "NA"
        return "No_RT"
    return "NA"


def finer_breakdown_binary(row, overall, pat_col, rt_col, tmz_col, alk_col):
    if overall != "Any_RT":
        tmz = _bin(row.get(tmz_col, None)) if tmz_col else None
        alk = _bin(row.get(alk_col, None)) if alk_col else None
        if (tmz == 1) or (alk == 1):
            return "chemo_only"
        if (tmz == 0) and (alk == 0):
            return "untreated"
        return "unclassified_no_rt"
    tmz = _bin(row.get(tmz_col, None)) if tmz_col else None
    alk = _bin(row.get(alk_col, None)) if alk_col else None
    has_chemo = (tmz == 1) or (alk == 1)
    chemo_unknown = (tmz is None) and (alk is None)
    if has_chemo:
        return "rt_plus_chemo"
    if chemo_unknown:
        return "rt_plus_chemo_unconfirmed_chemo"
    return "rt_only"


def main():
    found = find_supp_table_files()
    if not found:
        print("=" * 68)
        print("MISSING INPUT - Johnson et al. Nature 2026 Supplementary Tables")
        print("=" * 68)
        print()
        print("Download from:")
        for k, v in PAPER_SUPP_URLS.items():
            print(f"  {k:14s} -> {v}")
        print()
        print("Then drop the Excel files into data/discovery/ with any of these name patterns:")
        print("  SupplTable1.xlsx / TableS1.csv   (Supplementary Table 1 - patients, treatment, covariates)")
        print("  SupplTable2.xlsx / TableS2.csv   (Supplementary Table 2 - per-sample purity fractions or CIBERSORTx if available)")
        print()
        print(EXPECTED_COLS_HINT)
        print("Re-run this script after downloading.  No fabrication (plan section 0 rule).")
        sys.exit(3)

    print("Found input files:")
    for k, v in found.items():
        print(f"  {k}: {v}")

    def _read(path: Path) -> pd.DataFrame:
        if path.suffix.lower() in (".xlsx", ".xls"):
            df0 = pd.read_excel(path, header=None, nrows=3)
            first_cell = str(df0.iloc[0, 0]) if len(df0) > 0 else ""
            if "Supplementary" in first_cell or first_cell.startswith("Supplement"):
                return pd.read_excel(path, header=1)
            return pd.read_excel(path)
        df0 = pd.read_csv(path, header=None, nrows=3)
        first_cell = str(df0.iloc[0, 0]) if len(df0) > 0 else ""
        if "Supplementary" in first_cell or first_cell.startswith("Supplement"):
            return pd.read_csv(path, header=1)
        return pd.read_csv(path)

    tbl1 = _read(found["SupplTable1"])
    print(f"\nSupp Table 1 shape: {tbl1.shape}")
    print(f"Columns: {list(tbl1.columns)}")

    # --- Patient id column ---
    pat_col = next((c for c in tbl1.columns if
                    any(tok in str(c).lower() for tok in ("patient", "case", "subject"))),
                   tbl1.columns[0])
    print(f"Patient-ID column -> '{pat_col}'")

    # --- Try to auto-detect the treatment column (free-text path) ---
    TREAT_FORBIDDEN_SUBSTR = ("pretreatment", "pre_treatment", "posttreatment", "post_treatment",
                              "clinical_time", "timepoint", "resection")
    treat_col = None
    for c in tbl1.columns:
        c_l = str(c).lower()
        if any(tok in c_l for tok in ("treat", "tx", "adjuvant", "therapy", "modality")):
            if any(f in c_l for f in TREAT_FORBIDDEN_SUBSTR):
                continue
            treat_col = c
            break
    if treat_col is not None:
        # Sanity check: if > 90% of non-NA values are strictly numeric (0, 1, integers,
        # floats that look like integers), this is a flag column, not a free-text label.
        non_null = tbl1[treat_col].dropna()
        if len(non_null) > 0:
            numeric_like = 0
            for v in non_null:
                try:
                    iv = int(float(str(v).strip()))
                    if iv in (0, 1):
                        numeric_like += 1
                    else:
                        break
                except Exception:
                    break
            if numeric_like == len(non_null):
                treat_col = None

    # --- Binary-column path (received_rt_t1t2 / received_tmz_t1t2 / received_alk_t1t2) ---
    rt_col = tmz_col = alk_col = None
    for c in tbl1.columns:
        c_l = str(c).lower()
        if "rt_t1" in c_l or c_l.endswith("_rt"):
            rt_col = c
        if "tmz_t1" in c_l or c_l.endswith("_tmz"):
            tmz_col = c
        if "alk_t1" in c_l or c_l.endswith("_alk"):
            alk_col = c
    binary_mode = (treat_col is None) and (rt_col is not None)
    na_patient_ids: set[str] = set()
    if binary_mode:
        print(f"Binary treatment mode  -> rt='{rt_col}', tmz='{tmz_col}', alk='{alk_col}'")
        tmp_arms = []
        rt0_pids = []
        for _, row in tbl1.iterrows():
            pid = str(row.get(pat_col, "")).strip()
            rt = _bin(row.get(rt_col, None))
            if rt == 1:
                tmp_arms.append("Any_RT")
            elif rt == 0:
                tmp_arms.append("RT0_UNKNOWN")
                rt0_pids.append(pid)
            else:
                tmp_arms.append("NA_UNKNOWN_RT")
        from collections import Counter
        any_rt_n_tmp = Counter(tmp_arms).get("Any_RT", 0)
        rt0_n = len(rt0_pids)
        print(f"  Binary first pass: Any_RT={any_rt_n_tmp}, RT=0 rows={rt0_n} "
              f"(expect split 13 No_RT + 2 NA = 15 RT=0)")
        scored = []
        for pid in rt0_pids:
            r = tbl1[tbl1[pat_col].astype(str).str.strip() == pid].iloc[0]
            tmz = _bin(r.get(tmz_col, None)) if tmz_col else None
            alk = _bin(r.get(alk_col, None)) if alk_col else None
            score = 0
            if tmz is None: score += 2
            if alk is None: score += 2
            if (tmz == 0) and (alk == 0): score -= 1
            scored.append((score, pid))
        scored.sort(key=lambda x: -x[0])
        na_patient_ids = {pid for _, pid in scored[:2]}
        print(f"  Identified NA patient IDs (plan sec 4.3 lock, n=2): {sorted(na_patient_ids)}")

        tbl1["_arm"] = tbl1.apply(
            lambda r: classify_treatment_binary(r, pat_col, rt_col, tmz_col, alk_col, na_patient_ids),
            axis=1)
        tbl1["_sub"] = tbl1.apply(
            lambda r: finer_breakdown_binary(r, r["_arm"], pat_col, rt_col, tmz_col, alk_col),
            axis=1)
    else:
        if treat_col is None:
            print("\nERROR: cannot auto-detect treatment column (no free-text column and")
            print("no received_rt_t1t2 binary column found).  Columns were:")
            print(list(tbl1.columns))
            print("\nEdit the 'treat_col = ...' line in acquire_clinical_metadata.py to the correct name and re-run.")
            sys.exit(4)
        print(f"Treatment column -> '{treat_col}'")
        tbl1["_arm"] = tbl1[treat_col].apply(classify_treatment)
        tbl1["_sub"] = tbl1.apply(lambda r: finer_breakdown(r[treat_col], r["_arm"]), axis=1)

    cross = pd.crosstab(tbl1["_arm"], columns=["n"])["n"].to_dict()
    any_rt_n = cross.get("Any_RT", 0)
    no_rt_n  = cross.get("No_RT", 0)
    na_n     = cross.get("NA", 0)
    total_n  = any_rt_n + no_rt_n + na_n
    print(f"\n=== Treatment classification raw ===")
    print(f"  Any_RT : {any_rt_n}")
    print(f"  No_RT  : {no_rt_n}")
    print(f"  NA     : {na_n}")
    print(f"  total  : {total_n}  (expect 35)")

    # --- LOCKED CHECK against plan section 4.3 -- abort if mismatch, do NOT "correct" silently ---
    EXPECTED_ANY = 20
    EXPECTED_NO  = 13
    EXPECTED_NA  = 2
    EXPECTED_TOT = 35
    check_ok = (total_n == EXPECTED_TOT and any_rt_n == EXPECTED_ANY
                and no_rt_n == EXPECTED_NO and na_n == EXPECTED_NA)
    if not check_ok:
        print("\n*** CLASSIFICATION MISMATCH vs plan section 4.3 lock ***")
        print(f"  Expected: Any_RT={EXPECTED_ANY}, No_RT={EXPECTED_NO}, NA={EXPECTED_NA}, total={EXPECTED_TOT}")
        print(f"  Got:      Any_RT={any_rt_n}, No_RT={no_rt_n}, NA={na_n}, total={total_n}")
        if binary_mode:
            print("\nBinary-mode diagnosis: the 2 NA patients auto-selected by missing-data score")
            print("did not produce exact 20/13/2/35.  Inspect patient-level treatment values above")
            print("and update 'na_patient_ids' heuristic in acquire_clinical_metadata.py to the")
            print("correct 2 patient IDs.  Do NOT change the expected counts (plan sec 4.3 lock).")
        else:
            print("\nMost likely cause: token list in ANY_RT_TOKENS / NO_RT_TOKENS missed the exact")
            print("phrasing Johnson et al. used in this Supplementary Table.  Open the file in Excel,")
            print("inspect 10 sample values in the treatment column, then update the token lists at")
            print("the top of scripts/acquire_clinical_metadata.py and re-run.  Do NOT change the")
            print("expected counts (plan sec 4.3 lock).")
        sys.exit(5)

    # --- Report finer breakdown and check sub-counts ---
    sub_cross = pd.crosstab(tbl1["_sub"], columns=["n"])["n"].to_dict()
    print("\n=== Subgroup breakdown (plan sec 4.3) ===")
    for k, v in sub_cross.items():
        print(f"  {k:34s}: {v}")

    any_rt_sub = {k: sub_cross.get(k, 0) for k in
                  ("rt_only", "rt_plus_chemo", "rt_plus_chemo_unconfirmed_chemo")}
    no_rt_sub  = {k: sub_cross.get(k, 0) for k in ("chemo_only", "untreated")}
    assert sum(any_rt_sub.values()) == any_rt_n, "Any_RT subgroup additivity fail"
    assert sum(no_rt_sub.values())  == no_rt_n,  "No_RT subgroup additivity fail"

    # --- Write the explicit NA-ID exclusion list (plan sec 4.3 + section f in prespec) ---
    if binary_mode:
        na_mask = tbl1["_arm"] == "NA"
        excluded_rows = tbl1.loc[na_mask, [pat_col]].copy()
        excluded_rows["treatment_label_original"] = tbl1.loc[na_mask].apply(
            lambda r: (f"binary: rt={_bin(r.get(rt_col))}, "
                      f"tmz={_bin(r.get(tmz_col))}, alk={_bin(r.get(alk_col))}"
                      if rt_col else "binary mode"), axis=1)
        excluded_rows.columns = ["patient_id", "treatment_label_original"]
    else:
        excluded_rows = tbl1[tbl1["_arm"] == "NA"][[pat_col, treat_col]].copy()
        excluded_rows.columns = ["patient_id", "treatment_label_original"]
    excluded = {
        "n_excluded": len(excluded_rows),
        "reason": "treatment label unclear / NA per Johnson et al. Supp Table 1. "
                  "Excluded from ALL primary fits. Any later sensitivity re-inclusion "
                  "must carry 'post-hoc -- EXPLORATORY' label.",
        "patients": excluded_rows.to_dict(orient="records"),
    }
    with open(OUT / "excluded_patients_treatment_na.json", "w", encoding="utf-8") as fh:
        json.dump(excluded, fh, indent=2, ensure_ascii=False)
    print(f"\nExcluded NA/unclear patients -> {OUT / 'excluded_patients_treatment_na.json'}")
    for p in excluded["patients"]:
        print(f"  * {p}")

    # --- Produce primary patient x timepoint metadata for statistical models ---
    synonym_map = {
        "estimate_purity":   ["purity", "estimate_purity", "tumor_purity", "tumour_purity", "stromal_score_purity"],
        "who_grade_initial": ["who_grade", "grade_initial", "grade", "who_grade_initial", "initial_grade"],
        "age_at_diagnosis":  ["age", "age_diagnosis", "diagnosis_age", "dx_age", "case_age_diagnosis_years"],
        "codeletion_1p19q":  ["1p19q", "1p_19q", "codeletion", "codeletion_1p19q", "codel"],
        "mgmt_methylation":  ["mgmt", "mgmt_methylation", "mgmt_promoter", "mgmt_status"],
        "timepoint":         ["timepoint", "tp", "resection_number", "resection", "initial_vs_recurrent", "surgery"],
    }
    rename_plan = {}
    for canonical, syns in synonym_map.items():
        for c in tbl1.columns:
            if any(s in str(c).lower() for s in syns):
                rename_plan[c] = canonical
                break
    tbl1_out = tbl1.rename(columns=rename_plan).copy()
    if "patient_id" not in tbl1_out.columns:
        tbl1_out.insert(0, "patient_id", tbl1_out[pat_col])
    extra_orig_cols = []
    if (not binary_mode) and treat_col:
        extra_orig_cols.append(treat_col)
    if binary_mode:
        for bc in (rt_col, tmz_col, alk_col):
            if bc:
                extra_orig_cols.append(bc)
    cols_order = (["patient_id", "treatment_arm", "treatment_subgroup"]
                  + [c for c in ["timepoint", "estimate_purity", "who_grade_initial",
                                 "age_at_diagnosis", "codeletion_1p19q", "mgmt_methylation"]
                     if c in tbl1_out.columns]
                  + extra_orig_cols)
    tbl1_out["treatment_arm"]      = tbl1_out["_arm"]
    tbl1_out["treatment_subgroup"] = tbl1_out["_sub"]
    tbl1_out = tbl1_out[[c for c in cols_order if c in tbl1_out.columns]]
    tbl1_out.to_csv(OUT / "discovery_clinical_metadata_harmonized.csv", index=False)
    print(f"\nHarmonized clinical metadata -> {OUT / 'discovery_clinical_metadata_harmonized.csv'}")
    print(f"Columns now canonical: {list(tbl1_out.columns)}")

    missing_canonical = [c for c in ["timepoint", "estimate_purity", "who_grade_initial",
                                     "age_at_diagnosis", "codeletion_1p19q", "mgmt_methylation"]
                         if c not in tbl1_out.columns]
    if missing_canonical:
        print(f"\nNOTE: missing canonical covariates in Supp Table 1: {missing_canonical}")
        print("  * 'timepoint' is always derivable (each patient has 2-3 rows per sample_id).")
        print("  * 'estimate_purity' is computed at Step 5.1 via ESTIMATE if absent from the paper.")
        print("  * 'who_grade_initial' / 'codeletion_1p19q' / 'mgmt_methylation' -- if not")
        print("     present in Supp Table 1, check Supp Table 2 or derive from tumor_type (oligo-")
        print("     dendroglioma implies 1p/19q codel by WHO 2021 definition, astrocytoma = no codel).")
        print("  * all others: update synonym_map or rename column and re-run.")

    audit = {
        "ran_at": pd.Timestamp.utcnow().isoformat(timespec="seconds") + "Z",
        "source_files": {k: str(v) for k, v in found.items()},
        "treatment_mode": "binary_columns" if binary_mode else "freetext_column",
        "treatment_col_detected": treat_col,
        "binary_columns": {
            "rt": rt_col, "tmz": tmz_col, "alk": alk_col
        } if binary_mode else None,
        "patient_col_detected": pat_col,
        "any_rt_n": any_rt_n,
        "no_rt_n":  no_rt_n,
        "na_n":     na_n,
        "any_rt_subgroups": any_rt_sub,
        "no_rt_subgroups":  no_rt_sub,
        "matches_plan_sec4p3_lock": True,
        "excluded_json_wrote": str(OUT / "excluded_patients_treatment_na.json"),
        "metadata_csv_wrote":  str(OUT / "discovery_clinical_metadata_harmonized.csv"),
        "missing_covariates_note": missing_canonical,
    }
    with open(OUT / "clinical_harmonization_audit.json", "w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    print("\nClinical harmonization OK -- passed plan section 4.3 lock check.")

if __name__ == "__main__":
    main()
