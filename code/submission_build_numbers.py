# -*- coding: utf-8 -*-
"""Stage 1 - single source of truth for every number used in the submission package.

READ-ONLY with respect to the project.  Never imports a pipeline script
(importing acquire_snrnaseq_discovery.py truncates logs/acquire_snrnaseq.log).

Derives every value from the existing artefacts in outputs/ and writes:
  submission/supplementary/results_numbers.json
  submission/supplementary/numbers_provenance.csv

Prints an expected-vs-derived difference table and STOPS (non-zero exit) if any
of the three pre-declared sanity sentinels disagree.
"""
import csv
import hashlib
import io
import json
import os
import sys
import time

ROOT = r"<PROJECT_ROOT>"
os.chdir(ROOT)
INT = os.path.join(ROOT, "outputs", "intermediate")
OUT = os.path.join(ROOT, "submission")
SUP = os.path.join(OUT, "supplementary")

NUM = {}          # dotted name -> value
PROV = []         # provenance rows


def put(name, value, src, key=""):
    NUM[name] = value
    PROV.append({"name": name, "value": value, "source_file": src, "source_key": key})
    return value


def jload(rel):
    with io.open(os.path.join(ROOT, rel), "r", encoding="utf-8") as f:
        return json.load(f)


def rd_csv(rel):
    with io.open(os.path.join(ROOT, rel), "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def num(x):
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() in ("nan", "none", "na"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- short names
PRIMARY_SHORT = ["CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
                 "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils"]
PRIMARY_LONG = {
    "CD8_T": "CD8_T_cells", "CD4_conv": "CD4_T_cells_conventional",
    "Treg_FOXP3": "T_regulatory_cells_FOXP3", "NK": "NK_cells",
    "B": "B_cells", "Plasma": "Plasma_cells",
    "Microglia_P2RY12_TMEM119": "Microglia_P2RY12_TMEM119",
    "MoMF_CD163": "Monocyte_derived_macrophages_CD163",
    "DC": "Dendritic_cells", "Neutrophils": "Neutrophils",
}
# proportion-file column name -> short cell label
PROP_COL2SHORT = {
    "B_cells": "B", "CD8_T_cells": "CD8_T",
    "Microglia_P2RY12_TMEM119": "Microglia_P2RY12_TMEM119",
    "Monocyte_derived_macrophages_CD163": "MoMF_CD163", "NK_cells": "NK",
}
NON_IMMUNE_WINNERS = {"Astrocytes_non_malignant", "Malignant_astrocytes_oligo",
                      "Endothelial_cells", "Oligodendrocytes_OPC",
                      "Pericytes_Fibroblasts"}


def main():
    diffs = []
    sentinels = {}

    # =================================================== cohorts / design
    lmm = jload(r"outputs\intermediate\lmm_summary.json")
    put("cohorts.discovery.n_samples", lmm["n_samples"],
        r"outputs\intermediate\lmm_summary.json", "n_samples")
    put("cohorts.discovery.n_patients", lmm["n_patients"],
        r"outputs\intermediate\lmm_summary.json", "n_patients")

    cca = jload(r"outputs\intermediate\clinical_canonical_adapter_audit.json")
    put("cohorts.discovery.n_any_rt_patients", cca["counts"]["n_any_rt_patients"],
        r"outputs\intermediate\clinical_canonical_adapter_audit.json", "counts.n_any_rt_patients")
    put("cohorts.discovery.n_no_rt_patients", cca["counts"]["n_no_rt_patients"],
        r"outputs\intermediate\clinical_canonical_adapter_audit.json", "counts.n_no_rt_patients")

    spa = jload(r"outputs\intermediate\sn_pipeline_audit.json")
    put("cohorts.discovery.samples_input_snrna", spa["n_samples_input"],
        r"outputs\intermediate\sn_pipeline_audit.json", "n_samples_input")
    put("cohorts.discovery.cells_final_postcap", spa["n_cells_postcap_total"],
        r"outputs\intermediate\sn_pipeline_audit.json", "n_cells_postcap_total")
    put("cohorts.discovery.cells_precap_total", spa["n_cells_after_qc_precap_total"],
        r"outputs\intermediate\sn_pipeline_audit.json", "n_cells_after_qc_precap_total")
    put("cohorts.discovery.n_samples_capped", spa["n_samples_capped"],
        r"outputs\intermediate\sn_pipeline_audit.json", "n_samples_capped")
    put("cohorts.discovery.cap_cells_per_sample", spa["max_cells_per_sample"],
        r"outputs\intermediate\sn_pipeline_audit.json", "max_cells_per_sample")
    put("cohorts.discovery.n_hvg", spa["final_shape"][1],
        r"outputs\intermediate\sn_pipeline_audit.json", "final_shape[1]")
    put("cohorts.discovery.annotation_distribution", spa["annotation_distribution"],
        r"outputs\intermediate\sn_pipeline_audit.json", "annotation_distribution")
    put("cohorts.discovery.n_excluded_patients", cca["counts"]["n_excluded_na_treatment"],
        r"outputs\intermediate\clinical_canonical_adapter_audit.json",
        "counts.n_excluded_na_treatment")
    qcd = jload(r"outputs\intermediate\_tmp_qc_survivor_diag.json")
    _gb = [v["genes_before_qc"] for v in qcd["load_stats"].values()
           if isinstance(v, dict) and "genes_before_qc" in v]
    put("cohorts.discovery.n_genes_full", max(_gb),
        r"outputs\intermediate\_tmp_qc_survivor_diag.json",
        "max(load_stats[*].genes_before_qc) over loaded samples")

    lmmg = jload(r"outputs\intermediate\lmm_summary_glass.json")
    put("cohorts.glass.n_samples", lmmg["n_samples"],
        r"outputs\intermediate\lmm_summary_glass.json", "n_samples")
    put("cohorts.glass.n_cases", lmmg["n_patients"],
        r"outputs\intermediate\lmm_summary_glass.json", "n_patients")
    gpa = jload(r"outputs\intermediate\glass_pipeline_audit.json")
    put("cohorts.glass.anyrt_pairs", gpa["arm_by_timepoint"]["primary"]["Any_RT"],
        r"outputs\intermediate\glass_pipeline_audit.json", "arm_by_timepoint.primary.Any_RT")
    put("cohorts.glass.nort_pairs", gpa["arm_by_timepoint"]["primary"]["No_RT"],
        r"outputs\intermediate\glass_pipeline_audit.json", "arm_by_timepoint.primary.No_RT")
    put("cohorts.glass.purity_identifiable", gpa["purity"]["identifiable"],
        r"outputs\intermediate\glass_pipeline_audit.json", "purity.identifiable")

    cga = jload(r"outputs\intermediate\cgga_pipeline_audit.json")
    put("cohorts.cgga.version", cga["cgga_version"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "cgga_version")
    put("cohorts.cgga.samples_clinical_total", cga["samples_total_clinical"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "samples_total_clinical")
    put("cohorts.cgga.idh_mut_recurrent_total", cga["samples_idh_mut_recurrent_total"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "samples_idh_mut_recurrent_total")
    put("cohorts.cgga.n_analysed", cga["samples_idh_mut_recurrent_treatment_nonna"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "samples_idh_mut_recurrent_treatment_nonna")
    put("cohorts.cgga.anyrt_samples", cga["treatment_arm_counts"]["Any_RT"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "treatment_arm_counts.Any_RT")
    put("cohorts.cgga.nort_samples", cga["treatment_arm_counts"]["No_RT"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "treatment_arm_counts.No_RT")
    put("cohorts.cgga.n_methods", cga["n_methods"],
        r"outputs\intermediate\cgga_pipeline_audit.json", "n_methods")

    # =================================================== evaluability (discovery)
    prop = rd_csv(r"outputs\intermediate\sn_proportions_per_sample.csv")
    n_prop_samples = len(prop)
    put("evaluability.n_samples_with_proportions", n_prop_samples,
        r"outputs\intermediate\sn_proportions_per_sample.csv", "n rows")
    cols = [c for c in prop[0].keys() if c != "sample_id"]
    present = {}
    for c in cols:
        k = sum(1 for r in prop if num(r[c]) is not None and abs(num(r[c])) > 0)
        present[c] = k
    put("evaluability.present_counts_by_column", present,
        r"outputs\intermediate\sn_proportions_per_sample.csv", "non-zero cell counts per column")

    # evaluability is judged on the 10 PRE-SPEC PRIMARY types only
    ev = {}
    for s in PRIMARY_SHORT:
        longn = PRIMARY_LONG[s]
        cnt = present.get(longn, 0)
        ev[s] = {"long_name": longn, "present_samples": cnt,
                 "evaluable": cnt >= 10, "has_proportion_column": longn in cols}
    absent = [ev[s]["long_name"] for s in PRIMARY_SHORT if not ev[s]["has_proportion_column"]]
    put("evaluability.by_primary_cell_type", ev,
        r"outputs\intermediate\sn_proportions_per_sample.csv", "derived")
    put("evaluability.absent_cell_types_no_column", absent,
        r"outputs\intermediate\sn_proportions_per_sample.csv", "derived")
    put("evaluability.threshold_samples", 10,
        r"outputs\intermediate\sn_proportions_per_sample.csv", "rule: present in >=10 of 75")
    put("evaluability.n_primary_types",
        len(jload(r"outputs\sec4_prespec.json")["a_cell_types_primary"]),
        r"outputs\sec4_prespec.json", "len(a_cell_types_primary)")
    n_eval = sum(1 for s in PRIMARY_SHORT if ev[s]["evaluable"])
    put("evaluability.n_primary_evaluable", n_eval,
        r"outputs\intermediate\sn_proportions_per_sample.csv", "derived >=10 samples")
    put("evaluability.evaluable_cell_types",
        sorted(s for s in PRIMARY_SHORT if ev[s]["evaluable"]),
        r"outputs\intermediate\sn_proportions_per_sample.csv", "derived")
    sentinels["evaluable==2"] = (n_eval == 2)

    # =================================================== discovery LMM results
    lmm_rows = {r["cell_family"]: r for r in rd_csv(r"outputs\intermediate\lmm_interaction_summary.csv")}
    perm_rows = {r["cell_family"]: r for r in rd_csv(r"outputs\intermediate\snrna_direct_perm_pvalues.csv")}
    disc = {}
    for s in PRIMARY_SHORT:
        r = lmm_rows.get(s, {})
        p = perm_rows.get(s, {})
        has = num(r.get("beta_interaction")) is not None
        disc[s] = {
            "long_name": PRIMARY_LONG[s],
            "has_estimate": has,
            "beta": num(r.get("beta_interaction")),
            "ci_low": num(r.get("ci_low_interaction")),
            "ci_high": num(r.get("ci_high_interaction")),
            "asym_p": num(r.get("asym_p_interaction")),
            "bh_q": num(r.get("asym_q_interaction_bh")),
            "perm_p": num(p.get("perm_p")),
            "evaluable": ev[s]["evaluable"],
        }
    put("discovery.lmm_by_cell", disc,
        r"outputs\intermediate\lmm_interaction_summary.csv + snrna_direct_perm_pvalues.csv", "joined")
    put("discovery.n_primary_with_estimate", sum(1 for v in disc.values() if v["has_estimate"]),
        r"outputs\intermediate\lmm_interaction_summary.csv", "finite beta_interaction count")
    put("discovery.n_primary_converged", lmm["n_converged"],
        r"outputs\intermediate\lmm_summary.json", "n_converged")

    # =================================================== power
    ptab = {r["cell_type"]: r for r in rd_csv(r"outputs\intermediate\post_data_power_table.csv")}
    pa = jload(r"outputs\intermediate\post_data_power_audit.json")
    put("power.nsim", pa["nsim"], r"outputs\intermediate\post_data_power_audit.json", "nsim")
    put("power.seed", pa["seed"], r"outputs\intermediate\post_data_power_audit.json", "seed")
    put("power.icc_set", pa["icc_set"], r"outputs\intermediate\post_data_power_audit.json", "icc_set")
    put("power.f2_anchors", pa["f2_anchors"], r"outputs\intermediate\post_data_power_audit.json", "f2_anchors")
    put("power.design_anyrt", pa["n_rt"], r"outputs\intermediate\post_data_power_audit.json", "n_rt")
    put("power.design_nort", pa["n_nort"], r"outputs\intermediate\post_data_power_audit.json", "n_nort")
    put("power.n_patients", pa["n_patients"], r"outputs\intermediate\post_data_power_audit.json", "n_patients")
    # f2-anchored power is a design property -> identical across rows; take first row
    ref_row = ptab["CD8_T"] if "CD8_T" in ptab else list(ptab.values())[0]
    by_icc = {}
    for icc in ("icc40", "icc60", "icc70"):
        by_icc[icc] = {
            "f2028": num(ref_row["power_f2028_" + icc]),
            "f2040": num(ref_row["power_f2040_" + icc]),
        }
    put("power.f2anchored_by_icc", by_icc,
        r"outputs\intermediate\post_data_power_table.csv", "power_f2028_/power_f2040_ columns")
    put("power.headline_f2028_icc60", by_icc["icc60"]["f2028"],
        r"outputs\intermediate\post_data_power_table.csv", "power_f2028_icc60")
    put("power.headline_f2040_icc60", by_icc["icc60"]["f2040"],
        r"outputs\intermediate\post_data_power_table.csv", "power_f2040_icc60")
    put("power.headline_label", ref_row["headline_label_icc60"],
        r"outputs\intermediate\post_data_power_table.csv", "headline_label_icc60")

    # =================================================== GLASS
    glmm_rows = {r["cell_family"]: r for r in rd_csv(r"outputs\intermediate\lmm_interaction_summary_glass.csv")}
    glass = {}
    for s in PRIMARY_SHORT:
        r = glmm_rows.get(s, {})
        has = num(r.get("beta_interaction")) is not None
        glass[s] = {
            "long_name": PRIMARY_LONG[s], "has_estimate": has,
            "beta": num(r.get("beta_interaction")),
            "ci_low": num(r.get("ci_low_interaction")),
            "ci_high": num(r.get("ci_high_interaction")),
            "asym_p": num(r.get("asym_p_interaction")),
            "bh_q": num(r.get("asym_q_interaction_bh")),
        }
    put("glass.lmm_by_cell", glass,
        r"outputs\intermediate\lmm_interaction_summary_glass.csv", "all rows")
    put("glass.n_primary_estimable", sum(1 for v in glass.values() if v["has_estimate"]),
        r"outputs\intermediate\lmm_interaction_summary_glass.csv", "finite beta count")
    gdc = jload(r"outputs\intermediate\glass_direction_consistency.json")
    q_sig = gdc["glass_significant_cells"]
    put("glass.cells_bh_q_lt_0.05", q_sig,
        r"outputs\intermediate\glass_direction_consistency.json", "glass_significant_cells")
    put("glass.n_cells_bh_q_lt_0.05", gdc["n_glass_bh_q_lt_0.05"],
        r"outputs\intermediate\glass_direction_consistency.json", "n_glass_bh_q_lt_0.05")
    put("glass.n_primary_estimable_both", gdc["n_primary_estimable_both"],
        r"outputs\intermediate\glass_direction_consistency.json", "n_primary_estimable_both")
    sentinels["glass_sig=={B,Microglia}"] = (
        sorted(q_sig) == ["B", "Microglia_P2RY12_TMEM119"])
    # descriptive sign flip for the two discovery-evaluable cells
    flips = {}
    for s in ("NK", "Microglia_P2RY12_TMEM119"):
        b_d, b_g = disc[s]["beta"], glass[s]["beta"]
        flips[s] = {"discovery_beta": b_d, "glass_beta": b_g,
                    "flip_sign": (b_d is not None and b_g is not None and (b_d > 0) != (b_g > 0))}
    put("glass.sign_flip_vs_discovery_evaluable", flips,
        r"outputs\intermediate\lmm_interaction_summary*.csv", "descriptive, no agreement count")

    # =================================================== CGGA
    cstat = rd_csv(r"outputs\intermediate\cgga_crosssectional_stats.csv")
    endo = {}
    for r in cstat:
        if r["cell_family"] == "Endothelial":
            endo[r["method"]] = {"delta_median": num(r["delta_median_anyrt_minus_nort"]),
                                 "mw_p": num(r["mannwhitney_p_2tail"])}
    put("cgga.endothelial_by_method", endo,
        r"outputs\intermediate\cgga_crosssectional_stats.csv", "cell_family==Endothelial")
    prim_q = [num(r["mannwhitney_q_primary_bh"]) for r in cstat
              if r["is_primary"] == "True" and num(r["mannwhitney_q_primary_bh"]) is not None]
    put("cgga.best_primary_mw_q", min(prim_q) if prim_q else None,
        r"outputs\intermediate\cgga_crosssectional_stats.csv", "min mannwhitney_q_primary_bh")
    prim_p = [num(r["mannwhitney_p_2tail"]) for r in cstat
              if r["is_primary"] == "True" and num(r["mannwhitney_p_2tail"]) is not None]
    put("cgga.best_primary_mw_p", min(prim_p) if prim_p else None,
        r"outputs\intermediate\cgga_crosssectional_stats.csv", "min mannwhitney_p_2tail")
    put("cgga.n_primary_cells_any_method_q_lt_0.05",
        sum(1 for r in cstat if r["is_primary"] == "True"
            and num(r["mannwhitney_q_primary_bh"]) is not None
            and num(r["mannwhitney_q_primary_bh"]) < 0.05),
        r"outputs\intermediate\cgga_crosssectional_stats.csv", "derived")

    # =================================================== sensitivity
    ea = jload(r"outputs\intermediate\immune_subclust_audit_immune_subclust.json")
    put("sensitivity.compartment_cells", ea["compartment"]["n_cells"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "compartment.n_cells")
    put("sensitivity.ptprc_present", ea["compartment"]["ptprc_present"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "compartment.ptprc_present")
    put("sensitivity.resolutions", ea["predeclared"]["resolutions"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "predeclared.resolutions")
    put("sensitivity.evaluable_counts", ea["evaluable_counts"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "evaluable_counts")
    put("sensitivity.gate_pass", ea["gate"]["pass"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "gate.pass")
    put("sensitivity.gate_rule", ea["gate"]["rule"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "gate.rule")
    put("sensitivity.decision", ea["gate"]["decision"],
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json", "gate.decision")
    sentinels["sensitivity_gate_failed"] = (ea["gate"]["pass"] is False)

    contam = {}
    for res, cl in ea["clusters_per_resolution"].items():
        tot = len(cl)
        non = sum(1 for c in cl if c.get("winner") in NON_IMMUNE_WINNERS)
        contam[str(res)] = {"non_immune_clusters": non, "total_clusters": tot}
    put("sensitivity.contamination_by_resolution", contam,
        r"outputs\intermediate\immune_subclust_audit_immune_subclust.json",
        "clusters_per_resolution: non-immune winner count")
    evcsv = rd_csv(r"outputs\intermediate\immune_subclust_evaluability_immune_subclust.csv")
    ev_by_res = {}
    for r in evcsv:
        ev_by_res.setdefault(r["resolution"], []).append(
            {"cell_type": r["cell_type"], "evaluable": r["evaluable"] == "True"})
    put("sensitivity.evaluable_by_resolution", ev_by_res,
        r"outputs\intermediate\immune_subclust_evaluability_immune_subclust.csv", "evaluable column")

    # =================================================== pre-spec / provenance
    pj = os.path.join(ROOT, "outputs", "sec4_prespec.json")
    pm = os.path.join(ROOT, "outputs", "sec4_prespec.md")
    put("prespec.version", jload(r"outputs\sec4_prespec.json")["version"],
        r"outputs\sec4_prespec.json", "version")
    put("prespec.locked_at", jload(r"outputs\sec4_prespec.json")["locked_at"],
        r"outputs\sec4_prespec.json", "locked_at")
    put("prespec.sha256_json", sha256(pj), r"outputs\sec4_prespec.json", "sha256 of file")
    put("prespec.sha256_md", sha256(pm), r"outputs\sec4_prespec.md", "sha256 of file")
    put("permutation.n_permutations",
        jload(r"outputs\sec4_prespec.json")["d_permutation"]["n_permutations"],
        r"outputs\sec4_prespec.json", "d_permutation.n_permutations")
    put("permutation.seed",
        jload(r"outputs\sec4_prespec.json")["d_permutation"]["seed"],
        r"outputs\sec4_prespec.json", "d_permutation.seed")

    # =================================================== compare vs expected
    EXPECT = {
        "cohorts.discovery.n_samples": 71,
        "cohorts.discovery.n_patients": 33,
        "cohorts.discovery.n_any_rt_patients": 20,
        "cohorts.discovery.n_no_rt_patients": 13,
        "evaluability.n_samples_with_proportions": 75,
        "evaluability.n_primary_evaluable": 2,
        "evaluability.present_counts_by_column.B_cells": 6,
        "evaluability.present_counts_by_column.Monocyte_derived_macrophages_CD163": 5,
        "evaluability.present_counts_by_column.CD8_T_cells": 1,
        "evaluability.present_counts_by_column.Microglia_P2RY12_TMEM119": 72,
        "evaluability.present_counts_by_column.NK_cells": 71,
        "discovery.n_primary_with_estimate": 5,
        "discovery.lmm_by_cell.NK.beta": -0.07344359319202912,
        "discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.beta": 0.07577177148028941,
        "discovery.lmm_by_cell.NK.bh_q": 0.45622582241196,
        "discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.bh_q": 0.45622582241196,
        "discovery.lmm_by_cell.NK.perm_p": 0.35826417358264173,
        "discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.perm_p": 0.1927807219278072,
        "power.nsim": 800,
        "power.headline_f2028_icc60": 0.61875,
        "power.headline_f2040_icc60": 0.78,
        "cohorts.glass.n_samples": 108,
        "cohorts.glass.n_cases": 54,
        "cohorts.glass.anyrt_pairs": 45,
        "cohorts.glass.nort_pairs": 9,
        "glass.n_primary_estimable": 10,
        "glass.n_cells_bh_q_lt_0.05": 2,
        "sensitivity.compartment_cells": 36696,
        "sensitivity.gate_pass": False,
        "sensitivity.evaluable_counts.0.5": 1,
        "sensitivity.evaluable_counts.1.0": 2,
        "sensitivity.evaluable_counts.1.5": 2,
        "sensitivity.contamination_by_resolution.0.5.non_immune_clusters": 9,
        "sensitivity.contamination_by_resolution.0.5.total_clusters": 37,
        "sensitivity.contamination_by_resolution.1.0.non_immune_clusters": 16,
        "sensitivity.contamination_by_resolution.1.0.total_clusters": 58,
        "sensitivity.contamination_by_resolution.1.5.non_immune_clusters": 19,
        "sensitivity.contamination_by_resolution.1.5.total_clusters": 66,
        "cohorts.discovery.n_samples_capped": 67,
        "cohorts.discovery.cells_precap_total": 437454,
        "cohorts.discovery.cells_final_postcap": 221016,
        "cohorts.cgga.n_analysed": 129,
        "cohorts.cgga.anyrt_samples": 106,
        "cohorts.cgga.nort_samples": 23,
    }

    def get(dotted):
        """Resolve a dotted name against the flat NUM registry.  Keys may
        themselves contain dots (e.g. resolution labels like '0.5'), so after
        matching the longest registered prefix the remainder is re-joined
        greedily against the nested dict."""
        parts = dotted.split(".")
        for i in range(len(parts), 0, -1):
            prefix = ".".join(parts[:i])
            if prefix in NUM:
                cur = NUM[prefix]
                rest = parts[i:]
                while rest:
                    for j in range(len(rest), 0, -1):
                        key = ".".join(rest[:j])
                        if isinstance(cur, dict) and key in cur:
                            cur, rest = cur[key], rest[j:]
                            break
                    else:
                        return "<absent>"
                return cur
        return "<absent>"

    print("=" * 78)
    print("STAGE 1 - expected vs derived")
    print("=" * 78)
    for k, exp in EXPECT.items():
        got = get(k)
        if isinstance(exp, float) and isinstance(got, float):
            ok = abs(exp - got) < 1e-9
        else:
            ok = exp == got
        if not ok:
            diffs.append((k, exp, got))
            print("DIFF  %-64s expected=%s derived=%s" % (k, exp, got))
    if not diffs:
        print("no differences: every expected anchor matches a derived value")
    print()

    print("-" * 78)
    print("Sentinels")
    print("-" * 78)
    for k, v in sentinels.items():
        print("%-28s %s" % (k, "OK" if v else "**FAIL**"))
    print()

    # =================================================== write outputs
    os.makedirs(SUP, exist_ok=True)
    doc = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "project": "DGNC",
        "framing": "pre-specified negative-discovery report",
        "differences_vs_expected": [{"name": a, "expected": b, "derived": c} for a, b, c in diffs],
        "sentinels": sentinels,
        "values": NUM,
    }
    with io.open(os.path.join(SUP, "results_numbers.json"), "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
    with io.open(os.path.join(SUP, "numbers_provenance.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "value", "source_file", "source_key"])
        w.writeheader()
        for r in PROV:
            rr = dict(r)
            if isinstance(rr["value"], (dict, list)):
                rr["value"] = json.dumps(rr["value"], ensure_ascii=False)
            w.writerow(rr)
    print("wrote", os.path.join(SUP, "results_numbers.json"),
          "(%d values)" % len(NUM))
    print("wrote", os.path.join(SUP, "numbers_provenance.csv"),
          "(%d rows)" % len(PROV))

    if not all(sentinels.values()):
        print("!! SENTINEL FAILURE - stopping per task rule")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
