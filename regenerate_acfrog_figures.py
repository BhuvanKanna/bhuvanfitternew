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
WORM_MAP_XLSX = "Supplementary Data 1 trunc 20250702.xlsx"

# --- OE phenotype status, read from the proposal's Figure 2A -----------------
# The mcOverexpression "Any" column of Fig 2A defines OE sensitivity (shaded =
# mcOE phenotype present). These worm-gene lists were extracted from the vector
# cells of Fig 2A in ACFROG~1.PDF (mcOE-Any = dark-grey fill). The three
# merged human/worm rows resolve to single worm orthologs: KCNJ6/KCNJ15 -> irk-2,
# RRP1/RRP1B -> rrp-1, USP25/USP16 -> K02C4.3 (all OE-sensitive); pad-1 = n.d.
# (not overexpressed) and is excluded. This replaces the earlier reliance on
# genes_of_interest.json, which lacked the OE-tolerant (no-phenotype) genes.
SENS_GENES = [  # 24 OE-sensitive Hsa21 orthologs (mcOE phenotype present)
    "chaf-2", "cle-1", "dip-2", "dnsn-1", "eva-1", "pat-3", "irk-2",
    "Y54E10A.11", "mrps-6", "ncam-1", "F43G9.12", "pdxk-1", "pfk-1.1",
    "rcan-1", "rrp-1", "nrd-1", "Y105E8A.1", "hlh-34", "sod-1", "D1037.1",
    "unc-26", "Y74C10AL.2", "trpp-10", "K02C4.3",
]
TOL_GENES = [   # 23 OE-tolerant Hsa21 orthologs (no mcOE phenotype)
    "adr-2", "atp-3", "B0024.15", "cbs-1", "cct-8", "D1086.9", "H39E23.3",
    "igcm-1", "mbk-1", "F38B6.4", "stc-1", "itsn-1", "zig-10", "mrpl-39",
    "mtq-2", "pes-4", "pad-2", "ikb-1", "rnt-1", "set-29", "sod-5",
    "ubc-14", "wdr-4",
]


def build_worm_groups(table_index):
    """
    Resolve the Fig-2A worm gene names to fourparam-table transcript IDs
    (``w{n}_{transcript}``) using the ``GeneName``/``transcript`` columns of the
    Supplementary Data 1 mapping. Returns (sensitive_ids, tolerant_ids). Genes
    absent from the dataset (hlh-34, sod-5, H39E23.3) simply drop out.
    """
    import re
    xl = pd.read_excel(WORM_MAP_XLSX).dropna(subset=["transcript"])
    xl["tid"] = ("w" + xl["wwww"].astype(int).astype(str) + "_"
                 + xl["transcript"].astype(str))
    xl["seqname"] = xl["transcript"].str.replace(r"\.\d+$", "", regex=True)
    tabset = set(table_index)
    by_name, by_seq = {}, {}
    for _, r in xl.iterrows():
        if r["tid"] not in tabset:
            continue
        by_name.setdefault(str(r["GeneName"]), []).append(r["tid"])
        by_seq.setdefault(str(r["seqname"]), []).append(r["tid"])

    def resolve(names):
        out = []
        for n in names:
            ids = (by_name.get(n) or by_seq.get(n)
                   or by_seq.get(re.sub(r"\.\d+$", "", n)))
            if ids:
                out += ids
        return out

    return resolve(SENS_GENES), resolve(TOL_GENES)


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


def curve_xy(values, n=400, floor=EXCLUDE_AT_OR_BELOW):
    """Fit one gene and return (x grid, fitted 4-param curve, fit dict).

    ``floor`` sets the low-expression exclusion cutoff (values <= floor are
    dropped). Pass ``-inf`` to disable the exclusion entirely.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    v = v[v > floor]
    bf = BhuvanFitter(v, "g")
    d = bf.fit("fourparam")
    x = np.linspace(bf.hist_edges[0], bf.hist_edges[-1], n)
    return x, bf.fourparam_function(x), bf, d


# ===========================================================================
# FIGURE 3 -- worm
# ===========================================================================
def figure3(exclude=True):
    """Regenerate Figure 3.

    exclude=True  (default): use the ``> -1`` low-expression exclusion baked into
                  the fits (the pipeline's normal filter) -> acfrog_figure3_worm.png
    exclude=False: use the non-excluded table + include the -1 floor when
                   re-histogramming, for comparison -> acfrog_figure3_worm_nofilter.png
    The fit_success / 0<TI<1 / n_obs>=30 validity filters still apply in both.
    """
    table_csv = WORM_TABLE if exclude else "worm_fourparam_table.csv"
    floor = EXCLUDE_AT_OR_BELOW if exclude else -np.inf
    out_png = ("acfrog_figure3_worm.png" if exclude
               else "acfrog_figure3_worm_nofilter.png")
    table = pd.read_csv(table_csv).set_index("gene")

    # OE phenotype groups come from Figure 2A (see build_worm_groups), NOT from
    # genes_of_interest.json -- the json lacked the OE-tolerant (no-phenotype)
    # genes that form the "absent" control, which is exactly what panel D needs.
    sens_ids, tol_ids = build_worm_groups(table.index)

    vt = valid(table)
    present = [i for i in sens_ids if i in vt.index]   # OE-sensitive (Fig 2A)
    absent = [i for i in tol_ids if i in vt.index]     # OE-tolerant  (Fig 2A)

    ti = vt["truncationindex"]
    ti_present = ti.loc[present]
    ti_absent = ti.loc[absent]
    ti_all = ti  # all valid transcripts

    p_mwu = mannwhitneyu(ti_present, ti_absent, alternative="greater").pvalue
    STAR = 0.30  # a gene is "ceiling-like" if truncationindex exceeds this

    # load the worm expression matrix (genes as columns after transpose)
    master = pd.read_csv(WORM_CSV).set_index("strain").T

    def draw_curves(ax, ids, color, star=False, title="", nmax=30):
        ids = list(ids)[:nmax]
        for gid in ids:
            if gid not in master.columns:
                continue
            try:
                x, y, _bf, d = curve_xy(master[gid].values, floor=floor)
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

    # ---- Panel A: overlaid Gaussian fits, Hsa21 orthologs -----------------
    gsA = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[:, 0],
                                           hspace=0.4)
    axA1 = fig.add_subplot(gsA[0])
    axA2 = fig.add_subplot(gsA[1])
    draw_curves(axA1, present, "crimson", star=True,
                title="Hsa21 orthologs -- OE phenotype PRESENT (Fig 2A)\n"
                      "(★ = truncation index > 0.3)")
    draw_curves(axA2, absent, "0.4",
                title="Hsa21 orthologs -- OE phenotype ABSENT (Fig 2A)")
    axA1.text(-0.16, 1.2, "A", transform=axA1.transAxes,
              fontsize=15, fontweight="bold")

    # ---- Panel B: truncation index box-and-whisker ------------------------
    axB = fig.add_subplot(outer[:, 1])
    groups = [ti_absent.values, ti_present.values]
    box_cols = ["0.5", "crimson"]
    bp = axB.boxplot(
        groups, positions=[0, 1], widths=0.55, showfliers=False,
        patch_artist=True, medianprops=dict(color="black", lw=2),
        whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"),
    )
    for patch, col in zip(bp["boxes"], box_cols):
        patch.set_facecolor(col)
        patch.set_alpha(0.35)
        patch.set_edgecolor(col)
    # overlay individual genes as jittered points
    for xpos, data, col in [(0, ti_absent, "0.35"), (1, ti_present, "crimson")]:
        jit = xpos + (RNG.random(len(data)) - 0.5) * 0.28
        axB.scatter(jit, data, s=14, color=col, alpha=0.7,
                    edgecolors="white", linewidths=0.3, zorder=3)
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["Absent\n(tolerant)", "Present\n(sensitive)"],
                        fontsize=8)
    axB.set_ylabel("Truncation index", fontsize=9)
    axB.set_xlim(-0.6, 1.6)
    axB.set_ylim(-0.02, 0.7)
    axB.tick_params(labelsize=7)
    sig = "n.s." if p_mwu >= 0.05 else "*"
    axB.set_title(f"Mann-Whitney p = {p_mwu:.2f} ({sig})", fontsize=8)
    # label the top present genes
    for gid in ti_present.sort_values(ascending=False).index[:4]:
        axB.annotate(gid.split("_")[-1], (1.32, ti_present.loc[gid]),
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
        v = v[v > floor]
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

    # ---- Panel D: probability-normalized truncation-index histograms ------
    # Matches the original panel's "p" axis (per-bin fraction, sums to 1),
    # not a density (which blows up with narrow bins).
    axD = fig.add_subplot(outer[1, 2])
    bins = np.linspace(0, 1, 26)

    def prob_step(data, color, label):
        w = np.ones(len(data)) / len(data)
        axD.hist(data, bins=bins, weights=w, histtype="step",
                 color=color, lw=1.6, label=f"{label} (n={len(data)})")

    prob_step(ti_all, "0.5", "All transcripts")
    prob_step(ti_absent, "blue", "OE pheno absent")
    prob_step(ti_present, "red", "OE pheno present")
    axD.set_xlabel("Truncation index", fontsize=9)
    axD.set_ylabel("p (fraction of genes)", fontsize=9)
    axD.legend(fontsize=7)
    axD.tick_params(labelsize=7)
    axD.text(0.5, 0.55, r"$\mathrm{Trunc}=\dfrac{h_{right}}{h_{max}}$",
             transform=axD.transAxes, fontsize=11)
    axD.text(-0.2, 1.02, "D", transform=axD.transAxes,
             fontsize=15, fontweight="bold")

    filt = ("with > -1 exclusion filter" if exclude
            else "WITHOUT the > -1 exclusion filter (comparison)")
    fig.suptitle("Figure 3 (regenerated from repo data; OE groups from Fig 2A; "
                 f"{filt}): expression-distribution truncation in wild C. elegans",
                 fontsize=11, y=0.99)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Fig 3 [{'filtered' if exclude else 'NO filter'}]: "
          f"OE-sensitive={len(present)} OE-tolerant={len(absent)} "
          f"all={len(ti_all)}  MWU(present>absent) p={p_mwu:.3f}  -> {out_png}")


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
    figure3(exclude=True)          # normal pipeline (> -1 filter)
    figure3(exclude=False)         # comparison: no exclusion filter
    figure4()
    print("done")
