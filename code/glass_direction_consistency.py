#!/usr/bin/env python3
"""
Part 6 — GLASS validation, Step 4: direction-consistency QC table.

Per plan_implementation_status.md PART 6 (lines 452-472):

    "direction-consistency QC table: per K=10 PRIMARY cell type,
     sign(beta_discovery) == sign(beta_glass) -> TRUE/FALSE,
     GLASS one-sided p (same direction) reported."

Inputs
------
  outputs/intermediate/lmm_interaction_summary.csv       (discovery, snrna_direct)
  outputs/intermediate/lmm_interaction_summary_glass.csv (GLASS, glass_quantiseq)

Outputs
-------
  outputs/intermediate/glass_direction_consistency.csv
  outputs/intermediate/glass_direction_consistency.json

Honest-scope note (do NOT silently drop)
----------------------------------------
The discovery snRNA-direct proportions file
(outputs/intermediate/sn_proportions_per_sample.csv) contains only 10 cell-type
columns: 5 of them are README-listed K=10 PRIMARY members (NK, B, Microglia,
MoMF, DC). CD8_T / CD4_conv / Treg_FOXP3 / Plasma / Neutrophils were never
annotated in the discovery snRNA run, so plan gate P2-5 ("full cohort MUST
cover K=10 PRIMARY family") previously FAILED and those betas are structurally
NOT ESTIMABLE in discovery. P6-3 (>= 7/10 same sign) is therefore only
evaluable on the intersection. Both denominators are reported:
  - K=10 declared family (plan literal reading)
  - K_intersection = PRIMARY members estimable in BOTH cohorts
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INTER = PROJECT_ROOT / "outputs" / "intermediate"

DISC_PATH = INTER / "lmm_interaction_summary.csv"
GLASS_PATH = INTER / "lmm_interaction_summary_glass.csv"
OUT_CSV = INTER / "glass_direction_consistency.csv"
OUT_JSON = INTER / "glass_direction_consistency.json"

PRIMARY_CELLS = [
    "CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
    "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils",
]

# P6-3 threshold, stated literally in the plan.
P6_3_MIN_AGREE = 7
P6_1_MIN_PAIRS = 40


def _one_sided_same_direction(p_two: float, beta_disc: float, beta_glass: float) -> float:
    """One-sided p for H1: glass beta has the SAME sign as the discovery beta.

    NOTE: this is a *post-hoc directional* p (the direction is taken from the
    discovery estimate, not pre-registered per cell type), so it is reported
    for descriptive QC only and must NOT be read as a confirmatory p-value.
    """
    if not np.isfinite(p_two) or not np.isfinite(beta_disc) or not np.isfinite(beta_glass):
        return float("nan")
    if np.sign(beta_disc) == 0 or np.sign(beta_glass) == 0:
        return float("nan")
    return p_two / 2.0 if np.sign(beta_disc) == np.sign(beta_glass) else 1.0 - p_two / 2.0


def main() -> int:
    for p in (DISC_PATH, GLASS_PATH):
        if not p.exists():
            print(f"[!!] Missing required input: {p}")
            return 2

    disc = pd.read_csv(DISC_PATH).set_index("cell_family")
    glass = pd.read_csv(GLASS_PATH).set_index("cell_family")

    rows = []
    for cell in PRIMARY_CELLS:
        b_d = float(disc.loc[cell, "beta_interaction"]) if cell in disc.index else float("nan")
        b_g = float(glass.loc[cell, "beta_interaction"]) if cell in glass.index else float("nan")
        p_g = float(glass.loc[cell, "asym_p_interaction"]) if cell in glass.index else float("nan")
        q_g = float(glass.loc[cell, "asym_q_interaction_bh"]) if cell in glass.index else float("nan")

        disc_est = np.isfinite(b_d)
        glass_est = np.isfinite(b_g)
        both = disc_est and glass_est
        agree = bool(np.sign(b_d) == np.sign(b_g)) if both else None

        rows.append({
            "cell_family": cell,
            "estimable_discovery": disc_est,
            "beta_discovery": b_d,
            "estimable_glass": glass_est,
            "beta_glass": b_g,
            "glass_p_two_sided": p_g,
            "glass_q_bh_primary": q_g,
            "both_estimable": both,
            "sign_agreement": agree,
            "glass_p_one_sided_same_direction": (
                _one_sided_same_direction(p_g, b_d, b_g) if both else float("nan")
            ),
            "glass_q_lt_0.05": bool(np.isfinite(q_g) and q_g < 0.05),
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)

    n_declared = len(PRIMARY_CELLS)
    both_mask = out["both_estimable"]
    n_intersection = int(both_mask.sum())
    n_agree_intersection = int((out.loc[both_mask, "sign_agreement"]).sum())
    n_agree_declared = int((out["sign_agreement"] == True).sum())  # noqa: E712

    n_estimable_glass = int(out["estimable_glass"].sum())
    n_sig_glass = int(out["glass_q_lt_0.05"].sum())

    p6_3_verdict = (
        "PASS" if n_agree_declared >= P6_3_MIN_AGREE
        else ("NOT_EVALUABLE" if n_intersection < n_declared else "FAIL")
    )

    summary = {
        "artifact": "glass_direction_consistency.csv",
        "purpose": "Part 6 step 4 — P6-3 direction consistency (descriptive QC)",
        "n_primary_declared": n_declared,
        "n_primary_estimable_discovery": int(out["estimable_discovery"].sum()),
        "n_primary_estimable_glass": n_estimable_glass,
        "n_primary_estimable_both": n_intersection,
        "n_sign_agreement_declared_denominator": n_agree_declared,
        "n_sign_agreement_intersection_denominator": n_agree_intersection,
        "p6_3_threshold": f">= {P6_3_MIN_AGREE}/{n_declared}",
        "p6_3_verdict": p6_3_verdict,
        "p6_3_note": (
            "Discovery snRNA-direct annotation covers only 5/10 PRIMARY cell types "
            "(plan gate P2-5 FAILED); CD8_T/CD4_conv/Treg_FOXP3/Plasma/Neutrophils "
            "are structurally NOT ESTIMABLE in discovery, so the plan-literal "
            ">= 7/10 test cannot be evaluated. Intersection-denominator agreement "
            f"= {n_agree_intersection}/{n_intersection} (random expectation "
            f"{n_intersection/2:.1f}); this is below chance and must be reported "
            "as such, not spun."
        ),
        "n_glass_bh_q_lt_0.05": n_sig_glass,
        "glass_significant_cells": out.loc[out["glass_q_lt_0.05"], "cell_family"].tolist(),
        "disclaimer": (
            "glass_p_one_sided_same_direction is a post-hoc directional reading of "
            "the GLASS two-sided p (direction taken from discovery). It is NOT "
            "confirmatory and is provided for QC transparency only."
        ),
    }
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[OK] Wrote {OUT_CSV}")
    print(f"[OK] Wrote {OUT_JSON}")
    print()
    print(out.to_string(index=False))
    print()
    print(f"P6-3 (plan literal, >= {P6_3_MIN_AGREE}/{n_declared}): {p6_3_verdict}")
    print(f"  PRIMARY estimable in discovery : {summary['n_primary_estimable_discovery']}/{n_declared}")
    print(f"  PRIMARY estimable in GLASS     : {n_estimable_glass}/{n_declared}")
    print(f"  PRIMARY estimable in BOTH      : {n_intersection}/{n_declared}")
    print(f"  sign agreement (intersection)  : {n_agree_intersection}/{n_intersection}")
    print(f"  GLASS BH-FDR q<0.05 cells      : {summary['glass_significant_cells']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
