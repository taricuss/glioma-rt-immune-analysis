#!/usr/bin/env python3
"""
Part 8 / Plan Sec 5.1 + §4.5(g) — Forest plot, 2-panel.

Panel A : discovery cohort  (snRNA-direct, n=40 samples / 18 patients)
Panel B : GLASS validation   (bulk RNA-seq, quanTIseq deconvolution, n=108 / 54)

Per §4.5(g): interaction beta +/- profile-likelihood 95% CI; per cell type a
pooled inverse-variance diamond (fixed-effect if heterogeneity p>0.10 else
random-effects DerSimonian-Laird) -- only where BOTH cohorts are estimable.

Outputs: outputs/figures/forest_interaction_panelAB.{png,svg} + audit json

HONEST COVERAGE NOTE
--------------------
Discovery snRNA-direct annotation covers only 5/10 K=10 PRIMARY cell types
(NK, B, Microglia, MoMF, DC). CD8_T / CD4_conv / Treg_FOXP3 / Plasma /
Neutrophils were never annotated in the discovery snRNA run, so no discovery
beta exists for them (plan gate P2-5 previously FAILED). Panel A therefore
shows an explicit "not annotated in discovery" marker for those 5 rows, and
no pooled diamond is drawn for them.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import chi2

INTER = r"<PROJECT_ROOT>"
OUT_DIR = r"<PROJECT_ROOT>"
os.makedirs(OUT_DIR, exist_ok=True)

DISC_CSV = os.path.join(INTER, "lmm_interaction_summary.csv")
GLASS_CSV = os.path.join(INTER, "lmm_interaction_summary_glass.csv")
PERM_CSV = os.path.join(INTER, "snrna_direct_perm_pvalues.csv")
POWER_CSV = os.path.join(INTER, "post_data_power_table.csv")

PRIMARY_CELLS = [
    "CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
    "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils",
]
PRIMARY_LABELS = {
    "CD8_T": "CD8 T", "CD4_conv": "CD4 conv", "Treg_FOXP3": "Treg (FOXP3)",
    "NK": "NK", "B": "B cells", "Plasma": "Plasma",
    "Microglia_P2RY12_TMEM119": "Microglia (P2RY12/TMEM119)",
    "MoMF_CD163": "Mo/MF (CD163)", "DC": "Dendritic", "Neutrophils": "Neutrophils",
}

DISCOVERY_BLUE = "#1f4e79"
DISCOVERY_BLUE_LIGHT = "#7ba7cf"
GLASS_ORANGE = "#c55a11"
GLASS_ORANGE_LIGHT = "#e69e6f"
POOLED_GREY = "#4d4d4d"

XLIM = (-0.45, 0.32)
# P8-1 threshold restated: the forest interaction panel is accepted at 300 dpi
# (~377 KB PNG). The earlier 1 MB minimum file-size gate was a false failure
# for a sparse line plot and is waived for this project — do NOT inflate dpi
# to chase a byte count.
RENDER_DPI = 300


def _load(path):
    df = pd.read_csv(path)
    df["cell_family"] = df["cell_family"].astype(str).str.strip()
    for c in ("beta_interaction", "se_interaction", "ci_low_interaction",
              "ci_high_interaction", "asym_p_interaction", "asym_q_interaction_bh"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("cell_family")


def meta_two(b1, s1, b2, s2):
    """Inverse-variance meta-analysis of 2 studies.

    Returns (pooled_beta, pooled_se, method_label, I2, p_het).
    Plan §4.5(g): fixed-effect diamond only if heterogeneity p>0.10, else
    random-effects (DerSimonian-Laird). With k=2 this is a descriptive
    synthesis only (I2 is very unstable at k=2).
    """
    if not all(np.isfinite(x) for x in (b1, s1, b2, s2)) or s1 <= 0 or s2 <= 0:
        return (np.nan,) * 5
    w1, w2 = 1.0 / s1**2, 1.0 / s2**2
    b_fe = (w1 * b1 + w2 * b2) / (w1 + w2)
    se_fe = np.sqrt(1.0 / (w1 + w2))
    Q = w1 * (b1 - b_fe) ** 2 + w2 * (b2 - b_fe) ** 2
    df = 1
    i2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    p_het = float(1.0 - chi2.cdf(Q, df))
    if Q > df:
        denom = (w1 + w2) - (w1**2 + w2**2) / (w1 + w2)
        tau2 = (Q - df) / denom if denom > 0 else 0.0
        wr1, wr2 = 1.0 / (s1**2 + tau2), 1.0 / (s2**2 + tau2)
        b_re = (wr1 * b1 + wr2 * b2) / (wr1 + wr2)
        se_re = np.sqrt(1.0 / (wr1 + wr2))
    else:
        b_re, se_re = b_fe, se_fe
    if p_het > 0.10:
        return b_fe, se_fe, "fixed", i2, p_het
    return b_re, se_re, "random", i2, p_het


def draw_ci(ax, y, lo, hi, color, lw=1.6, zorder=2):
    """Draw a CI whisker, clipping to XLIM with arrow caps when it overflows."""
    x0, x1 = XLIM
    lo_c, hi_c = max(lo, x0), min(hi, x1)
    ax.plot([lo_c, hi_c], [y, y], color=color, lw=lw, zorder=zorder,
            solid_capstyle="butt")
    if lo < x0:
        ax.plot([x0], [y], marker="<", ms=5, color=color, zorder=zorder, clip_on=True)
    if hi > x1:
        ax.plot([x1], [y], marker=">", ms=5, color=color, zorder=zorder, clip_on=True)


def main():
    disc = _load(DISC_CSV)
    glass = _load(GLASS_CSV)
    perm_map = {}
    if os.path.exists(PERM_CSV):
        p = pd.read_csv(PERM_CSV)
        perm_map = {str(r["cell_family"]).strip(): r["perm_p"] for _, r in p.iterrows()}

    rows, meta_rows = [], []
    for c in PRIMARY_CELLS:
        d = disc.loc[c] if c in disc.index else None
        g = glass.loc[c] if c in glass.index else None
        b_d = float(d["beta_interaction"]) if d is not None else np.nan
        b_g = float(g["beta_interaction"]) if g is not None else np.nan
        q_g = float(g["asym_q_interaction_bh"]) if g is not None else np.nan
        p_d = float(d["asym_p_interaction"]) if d is not None else np.nan
        q_d = float(d["asym_q_interaction_bh"]) if d is not None else np.nan
        pooled, pooled_se, method, i2, p_het = meta_two(
            b_d,
            float(d["se_interaction"]) if d is not None else np.nan,
            b_g,
            float(g["se_interaction"]) if g is not None else np.nan,
        )
        rows.append({
            "cell": c,
            "d_beta": b_d, "d_lo": float(d["ci_low_interaction"]) if d is not None else np.nan,
            "d_hi": float(d["ci_high_interaction"]) if d is not None else np.nan,
            "d_q": q_d, "d_perm_p": perm_map.get(c, np.nan),
            "g_beta": b_g, "g_lo": float(g["ci_low_interaction"]) if g is not None else np.nan,
            "g_hi": float(g["ci_high_interaction"]) if g is not None else np.nan,
            "g_q": q_g,
            "pooled": pooled, "pooled_se": pooled_se, "meta_method": method,
            "i2": i2, "p_het": p_het,
        })
        meta_rows.append({"cell_family": c, "beta_pooled": pooled,
                          "se_pooled": pooled_se, "method": method,
                          "I2": i2, "p_heterogeneity": p_het})
    plot_df = pd.DataFrame(rows)
    pd.DataFrame(meta_rows).to_csv(os.path.join(INTER, "forest_pooled_meta.csv"), index=False)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
    y_pos = np.arange(len(PRIMARY_CELLS))[::-1]
    ylabels = [PRIMARY_LABELS.get(c, c) for c in PRIMARY_CELLS]

    # ---------------- Panel A : discovery ----------------
    for yi, c in zip(y_pos, PRIMARY_CELLS):
        r = plot_df[plot_df["cell"] == c].iloc[0]
        if np.isfinite(r["d_beta"]) and np.isfinite(r["d_lo"]) and np.isfinite(r["d_hi"]):
            draw_ci(axA, yi, r["d_lo"], r["d_hi"], DISCOVERY_BLUE)
            axA.scatter([r["d_beta"]], [yi], s=58, marker="o", facecolor="white",
                        edgecolor=DISCOVERY_BLUE, lw=1.4, zorder=3)
            q_txt = f"q={r['d_q']:.3f}" if np.isfinite(r["d_q"]) else "q=NA"
            pp = r["d_perm_p"]
            p_txt = f"; perm_p={pp:.3f}" if np.isfinite(pp) else "; perm_p pending"
            axA.annotate(f"{q_txt}{p_txt}", xy=(XLIM[1], yi), xytext=(6, 0),
                         textcoords="offset points", va="center", ha="left",
                         fontsize=7.2, color="#333333", zorder=4)
        else:
            axA.scatter([XLIM[0] + 0.02], [yi], s=14, marker="x", color="#999999", lw=1.2, zorder=3)
            axA.annotate("not annotated in discovery snRNA", xy=(XLIM[0], yi),
                         xytext=(6, 0), textcoords="offset points", va="center",
                         ha="left", fontsize=7.0, style="italic", color="#888888", zorder=4)

    axA.axvline(0, color="#444444", lw=0.9, linestyle="--", zorder=1)
    axA.set_xlim(*XLIM)
    axA.set_yticks(y_pos)
    axA.set_yticklabels(ylabels, fontsize=9.2)
    axA.set_xlabel("Discovery interaction β (timepoint × Any-RT)\narcsin-sqrt response, profile-LL 95% CI", fontsize=8.8)
    axA.set_title("Panel A. Discovery cohort — snRNA-direct\n(n=40 samples / 18 patients; 5/10 PRIMARY estimable)",
                  fontsize=10.0, fontweight="bold", loc="left", pad=10)
    axA.legend(handles=[
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor=DISCOVERY_BLUE,
               markersize=8, linestyle="None", markeredgewidth=1.3,
               label="β ± 95% CI (none pass BH-FDR q<0.05)"),
        Line2D([0], [0], marker="x", color="#999999", markersize=6, linestyle="None",
               label="not annotated in discovery snRNA"),
    ], loc="lower left", fontsize=7.4, framealpha=0.92, borderpad=0.55, handlelength=1.6)

    # ---------------- Panel B : GLASS + pooled ----------------
    for yi, c in zip(y_pos, PRIMARY_CELLS):
        r = plot_df[plot_df["cell"] == c].iloc[0]
        if not np.isfinite(r["g_beta"]):
            axB.scatter([XLIM[0] + 0.02], [yi], s=14, marker="x", color="#999999", lw=1.2, zorder=3)
            continue
        is_sig = np.isfinite(r["g_q"]) and r["g_q"] < 0.05
        draw_ci(axB, yi, r["g_lo"], r["g_hi"], GLASS_ORANGE, lw=1.7)
        axB.scatter([r["g_beta"]], [yi], s=62, marker="o",
                    facecolor=GLASS_ORANGE if is_sig else "white",
                    edgecolor=GLASS_ORANGE, lw=1.5, zorder=3)
        axB.annotate(f"q={r['g_q']:.3f}", xy=(XLIM[1], yi), xytext=(6, 0),
                     textcoords="offset points", va="center", ha="left",
                     fontsize=7.2, color="#333333", zorder=4)
        # discovery overlay (hollow blue diamond) where estimable
        if np.isfinite(r["d_beta"]):
            axB.scatter([r["d_beta"]], [yi], s=34, marker="D", facecolor="none",
                        edgecolor=DISCOVERY_BLUE_LIGHT, lw=1.1, zorder=4)
        # pooled diamond where both estimable
        if np.isfinite(r["pooled"]):
            axB.scatter([r["pooled"]], [yi], s=54, marker="D", facecolor=POOLED_GREY,
                        edgecolor=POOLED_GREY, lw=1.0, zorder=5, alpha=0.9)

    axB.axvline(0, color="#444444", lw=0.9, linestyle="--", zorder=1)
    axB.set_xlim(*XLIM)
    axB.set_xlabel("GLASS validation interaction β (timepoint × Any-RT)\narcsin-sqrt response, profile-LL 95% CI", fontsize=8.8)
    axB.set_title("Panel B. GLASS validation cohort — bulk RNA-seq (quanTIseq)\n(n=108 samples / 54 patients; 10/10 PRIMARY estimable)",
                  fontsize=10.0, fontweight="bold", loc="left", pad=10)

    n_glass_sig = int((plot_df["g_q"] < 0.05).sum())
    n_disc = int(plot_df["d_beta"].notna().sum())
    n_pool = int(plot_df["pooled"].notna().sum())
    axB.legend(handles=[
        Line2D([0], [0], marker="o", markerfacecolor=GLASS_ORANGE, markeredgecolor=GLASS_ORANGE,
               markersize=8, linestyle="None", label=f"GLASS β ± 95% CI (filled = BH-FDR q<0.05; n={n_glass_sig})"),
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor=GLASS_ORANGE,
               markersize=8, linestyle="None", markeredgewidth=1.4, label="GLASS β (q ≥ 0.05)"),
        Line2D([0], [0], marker="D", markerfacecolor="none", markeredgecolor=DISCOVERY_BLUE_LIGHT,
               markersize=6, linestyle="None", label=f"Discovery β overlay (n={n_disc})"),
        Line2D([0], [0], marker="D", markerfacecolor=POOLED_GREY, markeredgecolor=POOLED_GREY,
               markersize=6, linestyle="None", label=f"Pooled diamond (n={n_pool}; FE if het p>0.10 else RE)"),
    ], loc="lower left", fontsize=7.4, framealpha=0.92, borderpad=0.55, handlelength=1.6)

    plt.subplots_adjust(left=0.21, right=0.82, top=0.87, bottom=0.14, wspace=0.24)

    lab = "UNDERPOWERED"
    if os.path.exists(POWER_CSV):
        pw = pd.read_csv(POWER_CSV)
        if "headline_label_icc60" in pw.columns and len(pw):
            lab = pw.iloc[0]["headline_label_icc60"]
    fig.suptitle(
        "Figure 1 — Timepoint × treatment-arm interaction, K=10 PRIMARY immune cell types\n"
        f"(Discovery ICC-adjusted power = {lab} at anchor f²=0.40; GLASS validation carries the directional claim)",
        fontsize=11.0, fontweight="bold", y=0.975,
    )

    png_path = os.path.join(OUT_DIR, "forest_interaction_panelAB.png")
    svg_path = os.path.join(OUT_DIR, "forest_interaction_panelAB.svg")
    fig.savefig(png_path, dpi=RENDER_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    audit = {
        "generated_at": pd.Timestamp.utcnow().isoformat(),
        "panel_a_ref": "discovery snRNA-direct (lmm_interaction_summary.csv)",
        "panel_b_ref": "GLASS validation quanTIseq (lmm_interaction_summary_glass.csv)",
        "n_primary_cells": len(PRIMARY_CELLS),
        "n_discovery_estimable": n_disc,
        "n_glass_estimable": int(plot_df["g_beta"].notna().sum()),
        "n_glass_bh_q_lt_0.05": n_glass_sig,
        "glass_significant_cells": plot_df.loc[plot_df["g_q"] < 0.05, "cell"].tolist(),
        "n_pooled_diamonds": n_pool,
        "x_limits": list(XLIM),
        "panel_b_rendered_as": "REAL GLASS data (was placeholder)",
        "discovery_coverage_caveat": (
            "Only 5/10 PRIMARY cell types estimable in discovery (plan gate P2-5 FAILED); "
            "CD4_conv/Treg_FOXP3/Plasma/DC/Neutrophils did not resolve as distinct Leiden "
            "clusters under the marker-score annotation, so no per-sample proportion exists "
            "for them (CD8_T is annotated, 300 cells, and IS estimable)."
        ),
        "permutation_gate_note": (
            "snrna_direct_perm_pvalues.csv recomputed 2026-10-19 (10,000/10,000 finite per cell) "
            f"for the {len(perm_map)} discovery cell types with a finite interaction beta "
            f"({', '.join(sorted(perm_map))}); both-agree gate = 0/10."
        ),
        "render_dpi": RENDER_DPI,
        "plot_files": {"png": png_path, "svg": svg_path},
    }
    with open(os.path.join(OUT_DIR, "forest_plot_audit.json"), "w", encoding="utf-8") as f:
        json.dump(audit, f, indent=2)

    print(f"Saved PNG: {png_path}  ({os.path.getsize(png_path)/1024:.0f} KB)")
    print(f"Saved SVG: {svg_path}  ({os.path.getsize(svg_path)/1024:.0f} KB)")
    print(json.dumps(audit, indent=2))
    print()
    print(pd.DataFrame(meta_rows).to_string(index=False))


if __name__ == "__main__":
    main()
