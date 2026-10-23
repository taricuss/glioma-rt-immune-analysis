# -*- coding: utf-8 -*-
"""Stage 2 - figures for the DGNC submission package.

Reads ONLY existing aggregate artefacts (plus the already-published numbers
registry from Stage 1).  No pipeline is run or imported.

Outputs to submission/figures/:
  Figure1..Figure4, FigureS1..FigureS3   x   png(600dpi) / tiff(600dpi LZW) / pdf / svg
  figure_legends.docx, figure_legends.txt
"""
from __future__ import annotations

import csv
import io
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np

ROOT = r"<PROJECT_ROOT>"
INT = os.path.join(ROOT, "outputs", "intermediate")
OUT = os.path.join(ROOT, "submission")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

MM2IN = 1.0 / 25.4
W2 = 180 * MM2IN          # two-column width
W1 = 85 * MM2IN           # single-column width

# Okabe-Ito colourblind-safe palette
C = {
    "blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
    "red": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9",
    "yellow": "#F0E442", "grey": "#777777", "black": "#000000",
}

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,       # embed as TrueType -> text stays text
    "ps.fonttype": 42,
    "svg.fonttype": "none",   # keep text as text in SVG
    "figure.dpi": 120,
    "savefig.bbox": "standard",
})

LEGENDS = []   # (id, title, legend text, alt text)


def jload(rel):
    with io.open(os.path.join(ROOT, rel), "r", encoding="utf-8") as f:
        return json.load(f)


def rd_csv(rel):
    with io.open(os.path.join(ROOT, rel), "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def _nest(flat):
    """The registry is stored with dotted flat keys; rebuild the tree so the
    figure code can address values as nested dicts."""
    out = {}
    for k, v in flat.items():
        cur = out
        parts = k.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    return out


NUM = _nest(jload(r"submission\supplementary\results_numbers.json")["values"])


def reg(fig_id, title, legend, alt):
    LEGENDS.append({"id": fig_id, "title": title, "legend": legend, "alt": alt})


OVERFLOW = []


def _text_overflow(fig):
    """Programmatic stand-in for eyeballing each figure: every text artist must
    lie inside the figure canvas (catches clipped axis labels / annotations)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    bad = []
    for ax in fig.axes:
        for t in list(ax.texts) + [ax.xaxis.label, ax.yaxis.label]:
            if t is None or not t.get_text():
                continue
            bb = t.get_window_extent(renderer=r)
            if (bb.x0 < fb.x0 - 1 or bb.y0 < fb.y0 - 1
                    or bb.x1 > fb.x1 + 1 or bb.y1 > fb.y1 + 1):
                bad.append(repr(t.get_text()[:48]))
    return bad


def save_all(fig, stem, width_mm, height_mm):
    """Write png(600dpi) + tiff(600dpi LZW) + pdf + svg for one figure."""
    base = os.path.join(FIG, stem)
    fig.set_size_inches(width_mm * MM2IN, height_mm * MM2IN)
    bad = _text_overflow(fig)
    if bad:
        OVERFLOW.append((stem, bad))
    fig.savefig(base + ".png", dpi=600)
    fig.savefig(base + ".tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".svg")
    return base


# --------------------------------------------------------------- Figure 1
def fig1():
    fig, ax = plt.subplots(figsize=(W2, 4.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec, fs=7.5, bold=False):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                               linewidth=0.8, zorder=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, zorder=3,
                fontweight="bold" if bold else "normal", wrap=True)

    def arrow(x1, y1, x2, y2, style="-|>"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                     mutation_scale=7, linewidth=0.7,
                                     color="#444444", zorder=1))

    n = NUM["cohorts"]["discovery"]
    box(2, 84, 30, 12,
        "Discovery snRNA-seq (Johnson 2026)\nGSE326221 / Synapse syn17038081\n"
        "%d samples loaded, 0 QC skips" % n["samples_input_snrna"],
        "#eaf3fb", C["blue"])
    box(36, 84, 28, 12,
        "Cell QC + 2,000 HVG\nper-sample cap %d cells\n%d of %d samples capped"
        % (n["cap_cells_per_sample"], n["n_samples_capped"], n["samples_input_snrna"]),
        "#eaf3fb", C["blue"])
    box(68, 84, 30, 12,
        "Cluster-level argmax\nmarker-score annotation\nleiden r=0.8",
        "#eaf3fb", C["blue"])
    arrow(32, 90, 36, 90)
    arrow(64, 90, 68, 90)

    box(2, 64, 30, 13,
        "Modelled design joined to\nclinical canonical covariates\n"
        "%d samples / %d patients" % (n["n_samples"], n["n_patients"]),
        "#ffffff", "#333333")
    box(36, 64, 28, 13,
        "Cell-type proportions\nper sample (75 x 10)\n"
        "%d cells retained" % n["cells_final_postcap"], "#ffffff", "#333333")
    box(68, 64, 30, 13,
        "Evaluability rule\npresent in >=10 of 75 samples\n"
        "%d of 10 primary types evaluable" % NUM["evaluability"]["n_primary_evaluable"],
        "#fdf3e3", C["orange"])
    arrow(90, 84, 90, 77)
    arrow(17, 84, 17, 77)

    box(2, 44, 30, 13,
        "Interaction mixed model\narcsin-sqrt(fraction) ~ timepoint x arm\n"
        "+ covariates + (1|patient)  |  BH-FDR K=10",
        "#ffffff", "#333333")
    box(36, 44, 28, 13,
        "Patient-level permutation null\n10,000 shuffles, seed %d"
        % NUM["power"]["seed"], "#ffffff", "#333333")
    box(68, 44, 30, 13,
        "Post-data power check\nNSIM=%d, ICC 0.40/0.60/0.70\n"
        "f2=0.28: %.1f%%  f2=0.40: %.1f%%"
        % (NUM["power"]["nsim"], NUM["power"]["headline_f2028_icc60"] * 100,
           NUM["power"]["headline_f2040_icc60"] * 100), "#ffffff", "#333333")
    arrow(90, 64, 90, 57)
    arrow(17, 64, 17, 57)
    arrow(50, 64, 50, 57)

    box(2, 26, 94, 12,
        "Discovery result (PRIMARY): only NK and Microglia are evaluable; both null "
        "(BH q = %.3f); power BORDERLINE."
        % NUM["discovery"]["lmm_by_cell"]["NK"]["bh_q"],
        "#f2f2f2", "#333333")

    g = NUM["cohorts"]["glass"]
    cg = NUM["cohorts"]["cgga"]
    box(2, 6, 30, 15,
        "GLASS RNA silver-set (EXPLORATORY)\n%d samples / %d cases\n"
        "%d Any_RT pairs, %d No_RT pairs\nstock deconvolution signatures"
        % (g["n_samples"], g["n_cases"], g["anyrt_pairs"], g["nort_pairs"]),
        "#f6ecf4", C["purple"])
    box(35, 6, 30, 15,
        "CGGA mRNAseq_693 (EXPLORATORY)\n%d IDH-mutant recurrent samples\n"
        "%d Any_RT / %d No_RT (cross-sectional)"
        % (cg["n_analysed"], cg["anyrt_samples"], cg["nort_samples"]),
        "#f6ecf4", C["purple"])
    box(68, 6, 30, 15,
        "Immune-subclustering SENSITIVITY\ngate FAILED (%d/%d/%d of 10 at r=0.5/1.0/1.5)\n"
        "PTPRC absent from HVG space"
        % (NUM["sensitivity"]["evaluable_counts"]["0.5"],
           NUM["sensitivity"]["evaluable_counts"]["1.0"],
           NUM["sensitivity"]["evaluable_counts"]["1.5"]), "#f2f2f2", "#333333")
    arrow(50, 38, 50, 26)
    arrow(50, 26, 17, 21)
    arrow(50, 26, 50, 21)
    arrow(50, 26, 83, 21)

    fig.tight_layout(pad=0.3)
    save_all(fig, "Figure1", 180, 122)
    plt.close(fig)
    reg("Figure 1", "Study flow",
        "Study flow for the DGNC pre-specified analysis. Discovery single-nucleus "
        "RNA-seq (Johnson et al. 2026; GSE326221, Synapse syn17038081) was processed "
        "to %d cells over 2,000 highly variable genes with a per-sample cap of 3,000 "
        "cells; cell-type proportions were derived by cluster-level argmax of marker "
        "scores at Leiden resolution 0.8. Proportions were joined to harmonised "
        "clinical covariates yielding a modelled design of %d samples from %d patients. "
        "Evaluability required a cell type to be present in at least 10 of 75 samples, "
        "which %d of the 10 pre-specified primary types met. The primary test was a "
        "patient-random-intercept mixed model of arcsin-square-root fraction with a "
        "timepoint x treatment-arm interaction, corrected by Benjamini-Hochberg over "
        "the 10 primary types, with a 10,000-shuffle patient-level permutation null "
        "and a post-data power check (NSIM=%d). GLASS and CGGA are exploratory "
        "external arms; the pre-declared immune-subclustering sensitivity analysis "
        "failed its gate."
        % (NUM["cohorts"]["discovery"]["cells_final_postcap"],
           NUM["cohorts"]["discovery"]["n_samples"], NUM["cohorts"]["discovery"]["n_patients"],
           NUM["evaluability"]["n_primary_evaluable"], NUM["power"]["nsim"]),
        "Flow diagram with three discovery processing steps (sample loading, QC plus "
        "HVG selection, cluster-level annotation), a modelled design of 71 samples "
        "from 33 patients, an evaluability step retaining 2 of 10 primary cell types, "
        "and the primary mixed model with permutation null and power check. Three "
        "lower boxes show the exploratory GLASS and CGGA arms and the sensitivity "
        "analysis whose gate failed.")


# --------------------------------------------------------------- Figure 2
def fig2():
    prim = NUM["evaluability"]["by_primary_cell_type"]
    order = ["CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
             "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils"]
    labels = {"CD8_T": "CD8 T", "CD4_conv": "CD4 conv", "Treg_FOXP3": "Treg",
              "NK": "NK", "B": "B cells", "Plasma": "Plasma",
              "Microglia_P2RY12_TMEM119": "Microglia",
              "MoMF_CD163": "Mo/MF (CD163)", "DC": "Dendritic",
              "Neutrophils": "Neutrophils"}
    thr = NUM["evaluability"]["threshold_samples"]
    vals = [prim[s]["present_samples"] for s in order]
    cols = [C["green"] if prim[s]["evaluable"] else C["grey"] for s in order]

    fig, axes = plt.subplots(1, 2, figsize=(W2, 3.4),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, color=cols, height=0.66, edgecolor="white", linewidth=0.4)
    ax.axvline(thr, color=C["red"], linestyle="--", linewidth=0.9)
    ax.text(thr + 1.5, y[0] + 0.35, "evaluability threshold = %d" % thr,
            color=C["red"], fontsize=7, va="bottom")
    for yy, v in zip(y, vals):
        ax.text(v + 1.2, yy, str(v), va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[s] for s in order])
    ax.set_xlabel("samples with a non-zero proportion (of 75)")
    ax.set_xlim(0, 80)
    ax.set_ylim(-0.7, len(order) - 0.2)
    ax.text(-0.16, 1.03, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # panel B - per-sample distribution (aggregate display only) for the 2 evaluable
    ax2 = axes[1]
    prop = {r["sample_id"]: r for r in rd_csv(
        r"outputs\intermediate\sn_proportions_per_sample.csv")}
    clin = rd_csv(r"outputs\intermediate\clinical_canonical.csv")
    groups = [("primary", "No_RT"), ("primary", "Any_RT"),
              ("recurrent", "No_RT"), ("recurrent", "Any_RT")]
    glabels = ["Prim\nNo_RT", "Prim\nAny_RT", "Rec\nNo_RT", "Rec\nAny_RT"]
    cells = [("NK_cells", "NK"), ("Microglia_P2RY12_TMEM119", "Microglia")]
    width = 0.34
    xpos = np.arange(len(groups))
    for k, (col, lab) in enumerate(cells):
        data, ns = [], []
        for tp, arm in groups:
            sid = [c["sample_id"] for c in clin if c["timepoint"] == tp
                   and c["treatment_arm"] == arm]
            v = [float(prop[s][col]) for s in sid if s in prop]
            data.append(v)
            ns.append(len(v))
        off = (k - 0.5) * width
        bp = ax2.boxplot(data, positions=xpos + off, widths=width * 0.85,
                         patch_artist=True, showfliers=False,
                         medianprops=dict(color="white", linewidth=0.8),
                         whiskerprops=dict(linewidth=0.6),
                         capprops=dict(linewidth=0.6))
        for b in bp["boxes"]:
            b.set(facecolor=C["blue"] if k == 0 else C["orange"],
                  edgecolor="#333333", linewidth=0.5, alpha=0.9)
    for j, n in enumerate(ns):
        ax2.text(j, -0.045, "n=%d" % n, ha="center", va="top", fontsize=6.5)
    ax2.set_xticks(xpos)
    ax2.set_xticklabels(glabels)
    ax2.set_ylabel("per-sample proportion (arcsin-sqrt not applied)")
    ax2.set_ylim(-0.09, 1.0)
    ax2.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=C["blue"], edgecolor="#333333"),
                        plt.Rectangle((0, 0), 1, 1, facecolor=C["orange"], edgecolor="#333333")],
               labels=["NK cells", "Microglia"], frameon=False, loc="upper right", ncol=2)
    ax2.text(-0.18, 1.03, "B", transform=ax2.transAxes, fontsize=10, fontweight="bold")
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    fig.tight_layout(pad=0.4)
    save_all(fig, "Figure2", 180, 86)
    plt.close(fig)
    reg("Figure 2", "Cell-type evaluability in the discovery cohort",
        "Panel A: number of the 75 discovery samples in which each pre-specified "
        "primary cell type had a non-zero proportion. The dashed line marks the "
        "pre-declared evaluability threshold of 10 samples. Only NK cells (71 samples) "
        "and Microglia (72 samples) reached the threshold (green); B cells (6), "
        "monocyte-derived macrophages (5) and CD8 T cells (1) did not, and CD4 "
        "conventional T cells, regulatory T cells, plasma cells, dendritic cells and "
        "neutrophils had no proportion column at all. Panel B: aggregate per-sample "
        "proportion distributions for the two evaluable cell types, grouped by "
        "timepoint and treatment arm (box = interquartile range, centre line = median, "
        "whiskers = 1.5 x IQR, outliers not drawn); group sizes are printed beneath "
        "each box. No patient-level values are shown.",
        "Left panel: horizontal bar chart of sample counts per cell type with a "
        "dashed threshold line at ten samples; two bars exceed it. Right panel: paired "
        "boxplots for NK cells and Microglia across four timepoint-by-treatment groups "
        "with sample sizes annotated.")


# --------------------------------------------------------------- Figure 3
def fig3():
    disc = NUM["discovery"]["lmm_by_cell"]
    glass = NUM["glass"]["lmm_by_cell"]
    order = ["CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
             "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils"]
    labels = {"CD8_T": "CD8 T", "CD4_conv": "CD4 conv", "Treg_FOXP3": "Treg",
              "NK": "NK", "B": "B cells", "Plasma": "Plasma",
              "Microglia_P2RY12_TMEM119": "Microglia",
              "MoMF_CD163": "Mo/MF (CD163)", "DC": "Dendritic",
              "Neutrophils": "Neutrophils"}
    xlim = (-0.30, 0.30)

    fig, axes = plt.subplots(1, 2, figsize=(W2, 3.6), sharex=True)
    ypos = np.arange(len(order))[::-1]

    # ---- panel A: discovery
    ax = axes[0]
    ax.axvline(0, color="#888888", linewidth=0.7)
    for yy, s in zip(ypos, order):
        d = disc[s]
        if d["evaluable"] and d["beta"] is not None:
            lo, hi = d["ci_low"], d["ci_high"]
            col = C["blue"]
            clip = (lo < xlim[0]) or (hi > xlim[1])
            ax.errorbar([d["beta"]], [yy], xerr=[[d["beta"] - max(lo, xlim[0])],
                                                 [min(hi, xlim[1]) - d["beta"]]],
                        fmt="o", color=col, ecolor=col, markersize=4,
                        elinewidth=1.0, capsize=2, zorder=3)
            if clip:
                ax.plot([xlim[0] + 0.004], [yy], marker="<", color=col, markersize=4)
        else:
            ax.text(xlim[0] + 0.01, yy, "not evaluable (no estimate)",
                    fontsize=6.5, color="#999999", va="center")
    ax.set_yticks(ypos)
    ax.set_yticklabels([labels[s] for s in order])
    ax.set_xlim(*xlim)
    ax.set_xlabel("interaction coefficient (arcsin-sqrt fraction)")
    ax.text(-0.02, 1.03, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ---- panel B: GLASS + discovery overlay
    ax2 = axes[1]
    ax2.axvline(0, color="#888888", linewidth=0.7)
    for yy, s in zip(ypos, order):
        g = glass[s]
        if g["beta"] is not None:
            q = g["bh_q"]
            sig = (q is not None and q < 0.05)
            face = C["purple"] if sig else "white"
            ax2.errorbar([g["beta"]], [yy], xerr=[[g["beta"] - g["ci_low"]],
                                                  [g["ci_high"] - g["beta"]]],
                         fmt="o", color=C["purple"], ecolor=C["purple"],
                         markersize=4.4, elinewidth=1.0, capsize=2,
                         markerfacecolor=face, markeredgecolor=C["purple"], zorder=3)
    # discovery overlay (small hollow markers) for the two evaluable types
    for s in order:
        d = disc[s]
        if d["evaluable"] and d["beta"] is not None:
            yy = ypos[order.index(s)]
            ax2.plot([d["beta"]], [yy], marker="o", markersize=3,
                     markerfacecolor="none", markeredgecolor=C["blue"],
                     markeredgewidth=0.8, zorder=4)
    ax2.set_yticks(ypos)
    ax2.set_yticklabels([])
    ax2.set_xlim(*xlim)
    ax2.set_xlabel("interaction coefficient (arcsin-sqrt fraction)")
    ax2.text(-0.02, 1.03, "B", transform=ax2.transAxes, fontsize=10, fontweight="bold")
    handles = [
        Line2D([], [], marker="o", linestyle="none", color=C["purple"],
               markerfacecolor=C["purple"], markersize=4.4, label="GLASS, BH q < 0.05"),
        Line2D([], [], marker="o", linestyle="none", color=C["purple"],
               markerfacecolor="white", markersize=4.4, label="GLASS, BH q >= 0.05"),
        Line2D([], [], marker="o", linestyle="none", color=C["blue"],
               markerfacecolor="none", markersize=3, label="discovery (evaluable)"),
    ]
    ax2.legend(handles=handles, frameon=False, loc="lower right", fontsize=6.5)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    fig.tight_layout(pad=0.4)
    save_all(fig, "Figure3", 180, 91)
    plt.close(fig)
    reg("Figure 3", "Interaction estimates: discovery and exploratory GLASS",
        "Panel A: discovery-cohort timepoint x treatment-arm interaction coefficients "
        "with profile-likelihood 95% confidence intervals, shown only for the two "
        "evaluable primary cell types (NK, Microglia). The remaining eight primary "
        "types had no estimable discovery effect and are annotated as text rows with "
        "no glyph. Panel B: the same interaction estimated in the GLASS external "
        "cohort (EXPLORATORY), for all ten primary cell types; filled markers denote "
        "Benjamini-Hochberg q < 0.05 and hollow markers q >= 0.05. Small hollow blue "
        "markers overlay the two discovery-evaluable estimates for visual comparison. "
        "No pooled or meta-analytic diamonds are drawn, because pooling opposite-signed "
        "estimates would be read as agreement between cohorts, which the data do not "
        "support; "
        "this is a deviation from pre-specification section 4.5(g). Because no pooled "
        "estimate is drawn, no between-cohort heterogeneity statistic is reported. "
        "Confidence intervals that extend beyond the plotted range are marked with a "
        "triangle. Both panels share the x axis.",
        "Two side-by-side forest plots. Left: discovery interaction estimates with "
        "confidence intervals for NK and Microglia only, with eight grey text rows "
        "marked not evaluable. Right: GLASS estimates for all ten cell types, filled "
        "markers for two cell types with q below 0.05, with two small hollow overlay "
        "markers from the discovery cohort.")


# --------------------------------------------------------------- Figure 4
def fig4():
    evcsv = rd_csv(r"outputs\intermediate\immune_subclust_evaluability_immune_subclust.csv")
    res = ["0.5", "1.0", "1.5"]
    order = ["CD8_T_cells", "CD4_T_cells_conventional", "T_regulatory_cells_FOXP3",
             "NK_cells", "B_cells", "Plasma_cells", "Microglia_P2RY12_TMEM119",
             "Monocyte_derived_macrophages_CD163", "Dendritic_cells", "Neutrophils"]
    lab = {"CD8_T_cells": "CD8 T", "CD4_T_cells_conventional": "CD4 conv",
           "T_regulatory_cells_FOXP3": "Treg", "NK_cells": "NK", "B_cells": "B cells",
           "Plasma_cells": "Plasma", "Microglia_P2RY12_TMEM119": "Microglia",
           "Monocyte_derived_macrophages_CD163": "Mo/MF (CD163)",
           "Dendritic_cells": "Dendritic", "Neutrophils": "Neutrophils"}
    grid = {}
    for r in evcsv:
        grid[(r["cell_type"], r["resolution"])] = r["evaluable"] == "True"

    fig, axes = plt.subplots(1, 2, figsize=(W2, 3.3),
                            gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    ypos = np.arange(len(order))[::-1]
    for j, r in enumerate(res):
        for yy, ct in zip(ypos, order):
            on = grid.get((ct, r), False)
            ax.add_patch(Rectangle((j - 0.42, yy - 0.34), 0.84, 0.68,
                                   facecolor=C["green"] if on else "#f0f0f0",
                                   edgecolor="#bbbbbb", linewidth=0.4))
            ax.text(j, yy, "yes" if on else "no", ha="center", va="center",
                    fontsize=6.5, color="white" if on else "#777777",
                    fontweight="bold" if on else "normal")
    ax.set_xticks(range(len(res)))
    ax.set_xticklabels(["r = %s" % r for r in res])
    ax.set_yticks(ypos)
    ax.set_yticklabels([lab[c] for c in order])
    ax.set_xlim(-0.6, len(res) - 0.4)
    ax.set_ylim(-0.7, len(order) - 0.2)
    ax.text(-0.22, 1.03, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax2 = axes[1]
    cont = NUM["sensitivity"]["contamination_by_resolution"]
    n_imm = [cont[r]["non_immune_clusters"] for r in res]
    tot = [cont[r]["total_clusters"] for r in res]
    frac = [n / t for n, t in zip(n_imm, tot)]
    xp = np.arange(len(res))
    ax2.bar(xp, frac, color=C["red"], width=0.55, edgecolor="white", linewidth=0.4)
    for x, n, t, f in zip(xp, n_imm, tot, frac):
        ax2.text(x, f + 0.008, "%d of %d\nclusters" % (n, t), ha="center",
                 va="bottom", fontsize=6.5)
    ax2.set_xticks(xp)
    ax2.set_xticklabels(["r = %s" % r for r in res])
    ax2.set_ylabel("fraction of clusters with a\nnon-immune argmax label")
    ax2.set_ylim(0, max(frac) * 1.35)
    ax2.text(-0.28, 1.03, "B", transform=ax2.transAxes, fontsize=10, fontweight="bold")
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)

    fig.tight_layout(pad=0.4)
    save_all(fig, "Figure4", 180, 84)
    plt.close(fig)
    reg("Figure 4", "Pre-declared immune-subclustering sensitivity analysis",
        "Panel A: for each of the ten primary cell types and each pre-declared "
        "Leiden resolution (0.5, 1.0, 1.5), whether the cell type was evaluable in the "
        "immune-compartment subclustering (resolved to a cluster and present in at "
        "least 10 samples). Cell counts were 1, 2 and 2 at resolutions 0.5, 1.0 and "
        "1.5 against a gate requiring at least 6 of 10 at resolution 1.0 and at least "
        "5 of 10 at both 0.5 and 1.5, so the pre-declared gate failed. Panel B: "
        "contamination of the immune compartment, expressed as the number and "
        "fraction of subclusters whose argmax marker-score label was a non-immune "
        "cell type, at each resolution. PTPRC was absent from the highly variable gene "
        "space, so cells outside the annotated immune clusters could not be recovered.",
        "Left: a 10 by 3 grid of green and grey cells summarising evaluable status per "
        "cell type per resolution. Right: bar chart of the fraction of subclusters with "
        "a non-immune argmax label at three resolutions.")


# --------------------------------------------------------------- Figure S1
def figS1():
    ptab = {r["cell_type"]: r for r in rd_csv(
        r"outputs\intermediate\post_data_power_table.csv")}
    anchors = NUM["power"]["f2_anchors"]
    iccs = [("icc40", 0.40, C["sky"]), ("icc60", 0.60, C["blue"]),
            ("icc70", 0.70, C["grey"])]
    ref = ptab["CD8_T"]

    fig, ax = plt.subplots(figsize=(W1, 2.8))
    xs = np.array([0.0] + list(anchors))
    for key, icc, col in iccs:
        ys = np.array([0.05] + [float(ref["power_f2028_" + key]),
                                float(ref["power_f2040_" + key])])
        ax.plot(xs, ys * 100, marker="o", markersize=3.5, color=col,
                linewidth=1.0, label="ICC = %.2f" % icc)
    ax.axhline(80, color=C["red"], linestyle="--", linewidth=0.8)
    ax.text(0.281, 81.5, "80%", color=C["red"], fontsize=6.5)
    est = [s for s in ptab if ptab[s]["has_observed_fit"] == "True"]
    for s in est:
        f2o = float(ptab[s]["f2_obs_icc60"])
        if 0 < f2o <= max(anchors):
            po = float(ptab[s]["power_f2obs_icc60"]) * 100
            ax.plot([f2o], [po], marker="x", color=C["orange"], markersize=4)
    ax.plot([], [], marker="x", color=C["orange"], linestyle="none", markersize=4,
            label="observed f2 (ICC 0.60)")
    ax.set_xlabel("assumed interaction effect size (Cohen f2)")
    ax.set_ylabel("simulated power (%)")
    ax.set_xlim(-0.01, max(anchors) + 0.02)
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, loc="upper left", fontsize=6.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout(pad=0.3)
    save_all(fig, "FigureS1", 85, 71)
    plt.close(fig)
    reg("Figure S1", "Post-data power as a function of effect size and ICC",
        "Simulated power (NSIM=%d, seed %d) for the discovery timepoint x treatment "
        "interaction as a function of assumed effect size (Cohen f2) at assumed "
        "intraclass correlations of 0.40, 0.60 and 0.70. Power at f2=0.28 and 0.40 was "
        "%.1f%% and %.1f%% at ICC 0.60. Orange crosses mark the observed f2 of each "
        "primary cell type with an estimable discovery fit, evaluated at ICC 0.60; all "
        "observed effects fall below f2=0.28, and the design is labelled %s."
        % (NUM["power"]["nsim"], NUM["power"]["seed"],
           NUM["power"]["headline_f2028_icc60"] * 100,
           NUM["power"]["headline_f2040_icc60"] * 100,
           NUM["power"]["headline_label"]),
        "Line chart of simulated power against assumed effect size for three ICC "
        "values, with a dashed 80 percent reference line and orange cross markers for "
        "the observed effect sizes of five cell types.")

    # ----------------------------------------------------------- Figure S2
    endo = NUM["cgga"]["endothelial_by_method"]
    methods = [m for m in ("quantiseq", "epic", "mcpcounter") if m in endo]
    fig, ax = plt.subplots(figsize=(W1, 2.8))
    xp = np.arange(len(methods))
    vals = [endo[m]["delta_median"] for m in methods]
    ax.bar(xp, vals, color=[C["blue"], C["green"], C["orange"]][:len(methods)],
           width=0.55, edgecolor="white", linewidth=0.4)
    for x, m in zip(xp, methods):
        p = endo[m]["mw_p"]
        ax.text(x, max(vals) * 0.06, "p = %.3f" % p, ha="center", va="bottom",
                fontsize=6.5)
    ax.axhline(0, color="#666666", linewidth=0.7)
    ax.set_xticks(xp)
    ax.set_xticklabels([m.upper() if m == "epic" else
                        ("quanTIseq" if m == "quantiseq" else "MCP-counter")
                        for m in methods])
    ax.set_ylabel("median difference in transformed\nendothelial abundance (Any_RT - No_RT)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout(pad=0.3)
    save_all(fig, "FigureS2", 85, 71)
    plt.close(fig)
    reg("Figure S2", "CGGA cross-sectional endothelial result (EXPLORATORY)",
        "Exploratory cross-sectional comparison of transformed endothelial abundance "
        "between the Any_RT (%d) and No_RT (%d) IDH-mutant recurrent CGGA samples "
        "(mRNAseq_693) for the three deconvolution methods. Bars show the difference "
        "in median transformed abundance; Mann-Whitney two-tailed p values are printed "
        "above each bar. Endothelial was the only cell type with a consistently "
        "direction-consistent result across all three methods; the best primary "
        "Benjamini-Hochberg q value across all primary types was %.3f, so no primary "
        "cell-type result in this arm is significant. MCP-counter returns abundance "
        "scores rather than fractions and is used descriptively only."
        % (NUM["cohorts"]["cgga"]["anyrt_samples"], NUM["cohorts"]["cgga"]["nort_samples"],
           NUM["cgga"]["best_primary_mw_q"]),
        "Bar chart of median endothelial abundance difference between irradiated and "
        "non-irradiated CGGA samples for three deconvolution methods, with p values "
        "annotated above the bars.")

    # ----------------------------------------------------------- Figure S3
    g = NUM["cohorts"]["glass"]
    fig, ax = plt.subplots(figsize=(W1, 2.8))
    xp = np.arange(2)
    w = 0.36
    ax.bar(xp - w / 2, [g["anyrt_pairs"]] * 2, width=w, color=C["purple"],
           edgecolor="white", linewidth=0.4, label="Any_RT")
    ax.bar(xp + w / 2, [g["nort_pairs"]] * 2, width=w, color=C["grey"],
           edgecolor="white", linewidth=0.4, label="No_RT")
    for x, v in zip(xp - w / 2, [g["anyrt_pairs"]] * 2):
        ax.text(x, v + 1, str(v), ha="center", fontsize=6.5)
    for x, v in zip(xp + w / 2, [g["nort_pairs"]] * 2):
        ax.text(x, v + 1, str(v), ha="center", fontsize=6.5)
    ax.set_xticks(xp)
    ax.set_xticklabels(["primary", "recurrent"])
    ax.set_ylabel("paired samples")
    ax.set_ylim(0, max(g["anyrt_pairs"] * 1.25, 10))
    ax.legend(frameon=False, loc="upper right", fontsize=6.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout(pad=0.3)
    save_all(fig, "FigureS3", 85, 71)
    plt.close(fig)
    reg("Figure S3", "GLASS cohort composition (EXPLORATORY)",
        "Composition of the GLASS external cohort: %d paired samples from %d cases, "
        "with %d Any_RT and %d No_RT pairs at each timepoint. Only aggregate counts are "
        "shown. The small No_RT stratum limits the precision of the exploratory "
        "interaction estimate in this cohort."
        % (g["n_samples"], g["n_cases"], g["anyrt_pairs"], g["nort_pairs"]),
        "Grouped bar chart of GLASS sample counts by timepoint and treatment arm.")


def write_legends():
    # txt
    with io.open(os.path.join(FIG, "figure_legends.txt"), "w", encoding="utf-8") as f:
        for L in LEGENDS:
            f.write("%s. %s\n\n%s\n\nAlt text: %s\n\n%s\n"
                    % (L["id"], L["title"], L["legend"], L["alt"], "-" * 72))
    # docx
    from docx import Document
    from docx.shared import Pt
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)
    h = doc.add_heading("Figure legends", level=1)
    for L in LEGENDS:
        doc.add_heading("%s  %s" % (L["id"], L["title"]), level=2)
        p = doc.add_paragraph(L["legend"])
        p.paragraph_format.space_after = Pt(6)
        pa = doc.add_paragraph()
        r = pa.add_run("Alt text: " + L["alt"])
        r.italic = True
    doc.save(os.path.join(FIG, "figure_legends.docx"))


def verify():
    from PIL import Image
    print("-" * 78)
    print("Figure verification (dpi / pixel size)")
    print("-" * 78)
    bad = []
    for stem in ["Figure1", "Figure2", "Figure3", "Figure4",
                 "FigureS1", "FigureS2", "FigureS3"]:
        for ext in ("png", "tiff", "pdf", "svg"):
            p = os.path.join(FIG, stem + "." + ext)
            if not os.path.exists(p):
                bad.append((stem, ext, "MISSING"))
                print("MISSING %s.%s" % (stem, ext))
                continue
        for ext in ("png", "tiff"):
            p = os.path.join(FIG, stem + "." + ext)
            im = Image.open(p)
            dpi = im.info.get("dpi", (0, 0))
            ok = dpi and min(dpi) >= 599.5
            print("%-9s %-5s %5dx%-5d dpi=%.0f/%.0f  %s"
                  % (stem, ext, im.size[0], im.size[1], dpi[0], dpi[1],
                     "OK" if ok else "**BAD**"))
            if not ok:
                bad.append((stem, ext, "dpi=%s" % (dpi,)))
            im.close()
    return bad


def main():
    fig1()
    fig2()
    fig3()
    fig4()
    figS1()
    write_legends()
    bad = verify()
    print()
    print("-" * 78)
    print("Text-overflow check (programmatic clipping stand-in)")
    print("-" * 78)
    if OVERFLOW:
        for stem, items in OVERFLOW:
            print("OVERFLOW %-9s %s" % (stem, items))
    else:
        print("no text artist falls outside the figure canvas")
    print()
    if bad:
        print("!! figure problems:", bad)
        return 2
    print("all figures written to", FIG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
