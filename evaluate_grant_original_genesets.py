#!/usr/bin/env python
"""Replicate the grant's OWN gene-set definitions (not this repo's genome-wide
POS/TOL lists) and test them with this repo's truncationindex, to isolate
whether the gap between the grant's claimed effect and evaluate_censoring_hypothesis.py's
null result comes from a different truncation-index formula, or from different genes.

acfrog.md's preliminary-results paragraph specifically says the human evidence
came from "most of 1,400 genes predicted to be triplication sensitive (Collins
et al., 2022)" vs "774 triplication-tolerant genes (negative control)" -- NOT
the 44-gene curated OMIM/G2P list (positive_genes.txt) used everywhere else in
this repo's classifier work. Collins pTriplo >= 0.94 (this repo's own "high"
cutoff, see CLAUDE.md) yields 1,559 genes here -- already in the right ballpark
(if anything MORE genes than the grant's "1,400", so this isn't an underpowered
comparison). The tolerant list (positiveANDnegativeControlGenes.csv, 678 unique
symbols) is the same type of pangenome duplication-tolerant list the grant
describes (774) -- close enough in count/kind that gene selection isn't the
likely explanation on the TOL side.

This script tests the grant's own categorical framing directly: "most of the
triplication-sensitive genes have a ceiling" (truncationindex > 0) vs "most of
the tolerant genes do not" (truncationindex == 0) -- using MY truncationindex,
the real ~1,400-1,559-gene pTriplo positive set, and the same tolerant list.

Run: python evaluate_grant_original_genesets.py
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, ks_2samp, mannwhitneyu

from regenerate_grant_figures import _cereb_table, load_tol, valid

PTRIPLO_THRESHOLDS = [0.94, 0.95]
NAMED_EXAMPLES = ["APP", "SNCA", "PCSK9"]  # the grant's own named example genes
FIG_OUT = "grant_original_genesets_eval.png"


def load_ptriplo():
    d = pd.read_csv("dosage_sensitivity_collins2022.tsv.gz", sep="\t")
    d.columns = [c.lstrip("#") for c in d.columns]
    return d[["gene", "pTriplo"]].dropna()


def symbol_ti(symbols, raw, sym2ensg, filtered=False):
    """One truncationindex value per symbol (canonical = highest n_obs ENSG)."""
    v = valid(raw, filtered)
    rows = []
    for s in symbols:
        ensgs = [e for e in sym2ensg.get(s, []) if e in v.index]
        if not ensgs:
            continue
        sub = v.loc[ensgs]
        best = sub.sort_values("n_obs", ascending=False).iloc[0]
        rows.append({"symbol": s, "truncationindex": best["truncationindex"],
                     "n_obs": best["n_obs"]})
    return pd.DataFrame(rows)


def main():
    raw, sym2ensg = _cereb_table()
    ptriplo = load_ptriplo()
    tol_syms = load_tol()

    tol_df = symbol_ti(tol_syms, raw, sym2ensg, filtered=False)
    print(f"TOL (duplication-tolerant control): {len(tol_syms)} symbols requested, "
          f"{len(tol_df)} matched to a valid cerebellum fit (grant's own count: 774)")

    for thr in PTRIPLO_THRESHOLDS:
        pos_syms = set(ptriplo.loc[ptriplo["pTriplo"] >= thr, "gene"])
        pos_df = symbol_ti(pos_syms, raw, sym2ensg, filtered=False)
        print(f"\n{'='*70}\npTriplo >= {thr}  ({len(pos_syms)} genes genome-wide; "
              f"grant's own count: ~1,400)\n{'='*70}")
        print(f"POS matched to a valid cerebellum fit: {len(pos_df)}")

        pos_ceiling = (pos_df["truncationindex"] > 0).mean()
        tol_ceiling = (tol_df["truncationindex"] > 0).mean()
        print(f"Fraction with a ceiling (TI > 0): POS={pos_ceiling:.3f}  TOL={tol_ceiling:.3f}  "
              f"(grant's claim: 'most' POS have one, 'most' TOL do not)")

        table = [[int((pos_df["truncationindex"] > 0).sum()),
                 int((pos_df["truncationindex"] == 0).sum())],
                [int((tol_df["truncationindex"] > 0).sum()),
                 int((tol_df["truncationindex"] == 0).sum())]]
        odds, p_fisher = fisher_exact(table, alternative="greater")
        print(f"Fisher's exact (POS more likely to have a ceiling than TOL): "
              f"odds={odds:.3f}  p={p_fisher:.4g}")

        mwu_p = mannwhitneyu(pos_df["truncationindex"], tol_df["truncationindex"],
                             alternative="greater").pvalue
        ks_stat, ks_p = ks_2samp(pos_df["truncationindex"], tol_df["truncationindex"])
        print(f"Median TI -- POS={pos_df['truncationindex'].median():.4f}  "
              f"TOL={tol_df['truncationindex'].median():.4f}")
        print(f"Mann-Whitney U (POS > TOL): p={mwu_p:.4g}")
        print(f"KS test (POS vs TOL distributions differ): D={ks_stat:.3f}  p={ks_p:.4g}")

        print(f"\nNamed example genes (grant explicitly cites these as showing a ceiling):")
        for g in NAMED_EXAMPLES:
            row = symbol_ti([g], raw, sym2ensg, filtered=False)
            if len(row):
                ti = row.iloc[0]["truncationindex"]
                print(f"  {g:8s} truncationindex={ti:.4f}  {'(has a ceiling)' if ti > 0 else '(UNCAPPED -- no ceiling)'}")
            else:
                print(f"  {g:8s} not found / no valid fit in the cerebellum table")

    # Figure: CDF comparison at the primary threshold (0.94, this repo's convention)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pos_syms = set(ptriplo.loc[ptriplo["pTriplo"] >= 0.94, "gene"])
    pos_df = symbol_ti(pos_syms, raw, sym2ensg, filtered=False)
    fig, ax = plt.subplots(figsize=(7, 6))
    for df, color, label in [(pos_df, "crimson", f"POS (pTriplo>=0.94, n={len(pos_df)})"),
                             (tol_df, "steelblue", f"TOL (n={len(tol_df)})")]:
        x = np.sort(df["truncationindex"].values)
        y = np.arange(1, len(x) + 1) / len(x)
        ax.plot(x, y, lw=2, color=color, label=label)
    ax.set(xlabel="truncationindex", ylabel="cumulative fraction",
          title="Grant's own gene definitions (Collins pTriplo>=0.94 vs dup-tolerant),\n"
                "tested with this repo's truncationindex")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"\nSaved figure -> {FIG_OUT}")


if __name__ == "__main__":
    main()
