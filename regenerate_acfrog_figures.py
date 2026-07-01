# -*- coding: utf-8 -*-
"""
regenerate_acfrog_figures.py

Regenerate Figures 3 and 4 of the Pierce-lab OE-threshold grant proposal
(ACFROG~1.PDF -> acfrog.md) *faithfully from the data in this repo*, using the
single-source-of-truth fitting library ``bhuvanfitter.BhuvanFitter``.

Important honesty note
----------------------
These are FAITHFUL reproductions from the current repo data + the current
bounded ``truncationindex`` metric. They intentionally do NOT reproduce the
strong "expression ceiling" separation drawn in the original proposal figures:

  * Worm (Fig 3): OE-sensitive ``mco`` genes have only a slightly higher
    truncation index than the genomic background (Mann-Whitney n.s.), not the
    p < 0.004 in the proposal.
  * Cerebellum (Fig 4): the classic triplication genes (APP, SOD1, SNCA, ...)
    are mid-range Gaussians with tapering right tails; their truncation indices
    are ~0, consistent with the repo's own finding that truncationindex is
    ~uncorrelated with pTriplo. (The proposal's Fig 4 also used a different
    expression scaling: x-axis 0-8 for 212 people vs. this file's log2 ~6-15
    for 266 samples.)

Outputs (repo root):
    acfrog_figure3_worm.png
    acfrog_figure4_cerebellum.png

Run:
    python regenerate_acfrog_figures.py
"""

import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.stats import mannwhitneyu

from bhuvanfitter import BhuvanFitter

RNG = np.random.default_rng(0)
MIN_OBS = 30
EXCLUDE_AT_OR_BELOW = -1.0

WORM_CSV = "worm.csv"
WORM_TABLE = "worm_fourparam_table_excluded_at_or_below_-1.csv"
CEREB_CSV = "cerebellumlog2.csv"
CEREB_TABLE = "cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv"
GOI_JSON = "genes_of_interest.json"


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------
def valid(table):
    """The notebooks' shared analysis filter."""
    m = (
        table["fit_success"]
        & table["truncationindex"].notna()
        & (table["truncationindex"] > 0)
        & (table["truncationindex"] < 1)
        & (table["n_obs"] >= MIN_OBS)
    )
    return table[m]


def curve_xy(values, n=400):
    """Fit one gene and return (x grid, fitted 4-param curve, fit dict)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    v = v[v > EXCLUDE_AT_OR_BELOW]
    bf = BhuvanFitter(v, "g")
    d = bf.fit("fourparam")
    x = np.linspace(bf.hist_edges[0], bf.hist_edges[-1], n)
    return x, bf.fourparam_function(x), bf, d


# ===========================================================================
# FIGURE 3 -- worm
# ===========================================================================
def figure3():
    table = pd.read_csv(WORM_TABLE).set_index("gene")
    goi = json.load(open(GOI_JSON))

    mco_ids, lof_ids = set(), set()
    for k in ("mco_dev", "mco_behavior"):
        for _g, ids in goi[k].items():
            mco_ids.update(ids)
    for k in ("lof_dev", "lof_behavior"):
        for _g, ids in goi[k].items():
            lof_ids.update(ids)

    vt = valid(table)
    present = [i for i in mco_ids if i in vt.index]          # OE-sensitive
    lof_only = [i for i in (lof_ids - mco_ids) if i in vt.index]
    background = [i for i in vt.index if i not in mco_ids and i not in lof_ids]

    ti = vt["truncationindex"]
    ti_present = ti.loc[present]
    ti_bg = ti.loc[background]
    ti_all = ti  # all valid transcripts

    p_mwu = mannwhitneyu(ti_present, ti_bg, alternative="greater").pvalue
    STAR = 0.30  # a gene is "ceiling-like" if truncationindex exceeds this

    # load the worm expression matrix (genes as columns after transpose)
    master = pd.read_csv(WORM_CSV).set_index("strain").T

    def draw_curves(ax, ids, color, star=False, title="", nmax=25):
        ids = list(ids)[:nmax]
        for gid in ids:
            if gid not in master.columns:
                continue
            try:
                x, y, _bf, d = curve_xy(master[gid].values)
            except Exception:
                continue
            tival = table.loc[gid, "truncationindex"]
            ax.plot(x, y, color=color, lw=0.9, alpha=0.55)
            if star and pd.notna(tival) and tival > STAR:
                # star at the right ceiling
                ax.plot(d["right"], _bf.fourparam_function(d["right"]),
                        "*", color="red", ms=11, zorder=5)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Expression (TPM)", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        ax.tick_params(labelsize=7)

    fig = plt.figure(figsize=(13, 8.5))
    outer = gridspec.GridSpec(2, 3, width_ratios=[1.35, 0.75, 1.15],
                              height_ratios=[1, 1], wspace=0.32, hspace=0.42)

    # ---- Panel A: 2x2 overlaid Gaussian fits ------------------------------
    gsA = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[:, 0],
                                           wspace=0.38, hspace=0.5)
    axA1 = fig.add_subplot(gsA[0, 0])
    axA2 = fig.add_subplot(gsA[0, 1])
    axA3 = fig.add_subplot(gsA[1, 0])
    axA4 = fig.add_subplot(gsA[1, 1])
    draw_curves(axA1, lof_only, "0.4", title="Hsa21 orthologs\nOE pheno absent")
    draw_curves(axA2, present, "crimson", star=True,
                title="Hsa21 orthologs\nOE pheno present")
    bg_lo = ti_bg.sort_values().index[:25]
    bg_hi = ti_bg.sort_values(ascending=False).index[:25]
    draw_curves(axA3, bg_lo, "0.4", title="Other transcripts\n(low index)")
    draw_curves(axA4, bg_hi, "crimson", star=True,
                title="Other transcripts\n(high index)")
    axA1.text(-0.28, 1.28, "A", transform=axA1.transAxes,
              fontsize=15, fontweight="bold")

    # ---- Panel B: truncation index strip plot -----------------------------
    axB = fig.add_subplot(outer[:, 1])
    for xpos, data, col in [(0, ti_bg, "0.4"), (1, ti_present, "crimson")]:
        jit = xpos + (RNG.random(len(data)) - 0.5) * 0.35
        axB.scatter(jit, data, s=10, color=col, alpha=0.35,
                    edgecolors="none")
        axB.plot([xpos - 0.25, xpos + 0.25], [data.median()] * 2,
                 color="black", lw=2)
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["Absent\n(background)", "Present\n(mco)"], fontsize=8)
    axB.set_ylabel("Truncation index", fontsize=9)
    axB.set_ylim(0, 1)
    axB.set_title(f"p = {p_mwu:.2f} (n.s.)", fontsize=9)
    # label the top present genes
    for gid in ti_present.sort_values(ascending=False).index[:4]:
        axB.annotate(gid.split("_")[-1], (1.18, ti_present.loc[gid]),
                     fontsize=6, va="center")
    axB.text(-0.55, 1.02, "B", transform=axB.transAxes,
             fontsize=15, fontweight="bold")

    # ---- Panel C: 4 example genes -----------------------------------------
    gsC = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[0, 2],
                                           wspace=0.4, hspace=0.6)
    examples = ti_present.sort_values(ascending=False).index[:4]
    for ax_pos, gid in zip([(0, 0), (0, 1), (1, 0), (1, 1)], examples):
        ax = fig.add_subplot(gsC[ax_pos])
        v = master[gid].values.astype(float)
        v = v[np.isfinite(v)]
        v = v[v > EXCLUDE_AT_OR_BELOW]
        bf = BhuvanFitter(v, gid)
        d = bf.fit("fourparam")
        ax.bar(bf.hist_centers, bf.hist_counts,
               width=np.diff(bf.hist_edges).mean(),
               color="0.7", alpha=0.7)
        xs = np.linspace(bf.hist_edges[0], bf.hist_edges[-1], 300)
        ax.plot(xs, bf.fourparam_function(xs), color="red", lw=1.5)
        ax.set_title(f"{gid.split('_')[-1]}\nTI={d['truncationindex']:.2f}",
                     fontsize=7)
        ax.tick_params(labelsize=6)
    fig.text(0.685, 0.94, "C", fontsize=15, fontweight="bold")

    # ---- Panel D: normalized truncation-index histograms ------------------
    axD = fig.add_subplot(outer[1, 2])
    bins = np.linspace(0, 1, 26)
    axD.hist(ti_all, bins=bins, density=True, histtype="step",
             color="0.5", lw=1.5, label="All transcripts")
    axD.hist(ti_bg, bins=bins, density=True, histtype="step",
             color="blue", lw=1.5, label="OE pheno absent")
    axD.hist(ti_present, bins=bins, density=True, histtype="step",
             color="red", lw=1.5, label="OE pheno present")
    axD.set_xlabel("Truncation index", fontsize=9)
    axD.set_ylabel("Density", fontsize=9)
    axD.legend(fontsize=7)
    axD.tick_params(labelsize=7)
    axD.text(0.55, 0.55, r"$\mathrm{Trunc}=\dfrac{h_{right}}{h_{max}}$",
             transform=axD.transAxes, fontsize=11)
    axD.text(-0.2, 1.02, "D", transform=axD.transAxes,
             fontsize=15, fontweight="bold")

    fig.suptitle("Figure 3 (regenerated from repo data): expression-distribution "
                 "truncation in wild C. elegans", fontsize=12, y=0.99)
    fig.savefig("acfrog_figure3_worm.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Fig 3: present={len(present)} lof_only={len(lof_only)} "
          f"background={len(background)}  MWU p={p_mwu:.3f}")


# ===========================================================================
# FIGURE 4 -- cerebellum (GTEx)
# ===========================================================================
HSA21 = ["APP", "SOD1", "DYRK1A", "RCAN1", "SYNJ1", "ITSN1",
         "DSCAM", "CBS", "DONSON", "PCP4"]
TRIP = ["SNCA", "PCSK9", "PMP22", "SCN2A", "SCN3A", "MECP2",
        "KCNQ2", "CACNA1C"]
DISEASE = {"APP": "Alzheimer's", "SNCA": "Parkinson's",
           "PCSK9": "Hypercholest.", "PMP22": "Charcot-Marie-Tooth 1A",
           "SCN2A": "Seizure", "MECP2": "MECP2 dup."}


def figure4():
    ids = pd.read_csv(CEREB_CSV, usecols=["Name", "Description"])
    sym2name = dict(zip(ids["Description"], ids["Name"]))
    table = pd.read_csv(CEREB_TABLE).set_index("gene")

    wanted_syms = HSA21 + TRIP
    wanted_names = {sym2name[s] for s in wanted_syms if s in sym2name}

    # read the big matrix, keep only the rows we need, then transpose
    full = pd.read_csv(CEREB_CSV)
    sub = full[full["Name"].isin(wanted_names)].drop(columns=["Description"])
    master = sub.set_index("Name").T   # samples x (selected genes)
    del full

    def panel(ax, syms, title):
        # rank by truncation index so we can star the most ceiling-like
        tis = {}
        for s in syms:
            n = sym2name.get(s)
            if n in table.index:
                tis[s] = table.loc[n, "truncationindex"]
        star_syms = set(sorted(tis, key=lambda s: -(tis.get(s) or 0))[:3])
        for s in syms:
            n = sym2name.get(s)
            if n is None or n not in master.columns:
                continue
            try:
                x, y, bf, d = curve_xy(master[n].values)
            except Exception:
                continue
            ceiling = s in star_syms
            col = "red" if ceiling else "0.45"
            ax.plot(x, y, color=col, lw=1.4 if ceiling else 1.0,
                    alpha=0.85 if ceiling else 0.6)
            # label at the peak
            xpk = d["x0"]
            ax.annotate(s, (xpk, bf.fourparam_function(xpk)),
                        fontsize=7, color=col, fontweight="bold",
                        ha="center", va="bottom")
            if ceiling:
                ax.plot(d["right"], bf.fourparam_function(d["right"]),
                        "*", color="red", ms=12, zorder=6)
            if s in DISEASE:
                ax.annotate(DISEASE[s], (xpk, bf.fourparam_function(xpk)),
                            xytext=(0, 14), textcoords="offset points",
                            fontsize=6, color="0.3", ha="center")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Expression (log2 TPM)", fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        ax.tick_params(labelsize=7)

    fig, (axA, axB) = plt.subplots(2, 1, figsize=(7.5, 8.5))
    panel(axA, HSA21, "A  Example Hsa21 genes")
    panel(axB, TRIP, "B  Example triplication-sensitive genes")
    for ax, lab in [(axA, "A"), (axB, "B")]:
        ax.text(-0.1, 1.04, lab, transform=ax.transAxes,
                fontsize=15, fontweight="bold")
    fig.suptitle("Figure 4 (regenerated from repo data): cerebellum expression "
                 "distributions\n(red = highest truncation index; ceilings are "
                 "weak in this dataset)", fontsize=10, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig("acfrog_figure4_cerebellum.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Fig 4: genes plotted =", len(master.columns))


if __name__ == "__main__":
    figure3()
    figure4()
    print("done")
