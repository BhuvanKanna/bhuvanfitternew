#!/usr/bin/env python
"""Follow-up to the Enrichr confound check: TOL is dominated by testis-restricted
genes that were never a fair "does this gene tolerate overexpression" comparison
to POS's brain-leaning profile. This builds a TISSUE-MATCHED TOL set -- each POS
gene matched to its k nearest TOL neighbors specifically on brain-vs-rest tissue
expression (log2fc_brain_vs_nonbrain), the same nearest-neighbor design already
used for expression level in evaluate_censoring_hypothesis.py, just applied to
tissue profile instead of raw mean -- then re-runs both diagnostics that showed a
brain-driven signal (Enrichr enrichment, and the tissue-aware classifier ablation)
to see what survives once this specific confound is actually controlled for.

Run: python tissue_matched_reanalysis.py
"""
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate
from sklearn.neighbors import NearestNeighbors

from ablation_dosage_features import LEVEL, PURE_SHAPE
from enrichr_tissue_enrichment import (LIBRARIES, get_enrichment, submit_list,
                                       top_brain_hit)
from regenerate_grant_figures import load_tol
from train_dosage_classifier import load_labelled, make_pipeline
from train_tissue_aware_classifier import TISSUE, load_tissue_features

K_NEIGHBORS = 5
POS_FILE = "positive_genes_compiled.txt"
MATCHED_TOL_OUT = "tol_tissue_matched.txt"
SUMMARY_OUT = "tissue_matched_enrichr_summary.csv"
FIG_OUT = "tissue_matched_comparison.png"


def build_tissue_matched_tol(k=K_NEIGHBORS):
    tissue = load_tissue_features()  # index=symbol, includes log2fc_brain_vs_nonbrain
    pos = sorted({l.strip() for l in open(POS_FILE) if l.strip()})
    tol = sorted(load_tol())

    pos_t = tissue.reindex(pos).dropna(subset=["log2fc_brain_vs_nonbrain"])
    tol_t = tissue.reindex(tol).dropna(subset=["log2fc_brain_vs_nonbrain"])
    print(f"POS with tissue data: {len(pos_t)}/{len(pos)}   "
         f"TOL with tissue data: {len(tol_t)}/{len(tol)}")

    nn = NearestNeighbors(n_neighbors=min(k, len(tol_t))).fit(tol_t[["log2fc_brain_vs_nonbrain"]])
    _, idx = nn.kneighbors(pos_t[["log2fc_brain_vs_nonbrain"]])
    matched_idx = sorted(set(idx.ravel()))
    matched_tol = tol_t.iloc[matched_idx]

    print(f"Matched TOL: {len(matched_tol)} genes (from {len(tol_t)} total, k={k})")
    print(f"Median log2fc_brain_vs_nonbrain -- POS: {pos_t['log2fc_brain_vs_nonbrain'].median():.3f}  "
         f"all-TOL: {tol_t['log2fc_brain_vs_nonbrain'].median():.3f}  "
         f"matched-TOL: {matched_tol['log2fc_brain_vs_nonbrain'].median():.3f}")

    with open(MATCHED_TOL_OUT, "w", newline="\n") as f:
        for g in matched_tol.index:
            f.write(g + "\n")

    return pos_t, tol_t, matched_tol


def rerun_enrichr(pos_symbols, matched_tol_symbols):
    print("\n== Re-running Enrichr: POS vs tissue-matched TOL ==")
    results = []
    for name, genes in [("POS", pos_symbols), ("TOL_tissue_matched", matched_tol_symbols)]:
        uid = submit_list(list(genes), name)
        time.sleep(1)
        for lib in LIBRARIES:
            df = get_enrichment(uid, lib)
            df.to_csv(f"enrichr_{name}_{lib}.csv", index=False)
            top = df.iloc[0]
            brain_hit = top_brain_hit(df, top_n=50)
            print(f"  {name:20s} {lib:24s} top: {top['term'][:40]:40s} adj_p={top['adj_pvalue']:.2e}"
                 + (f"  | brain (rank {int(brain_hit['rank'])}): {brain_hit['term'][:35]} "
                    f"adj_p={brain_hit['adj_pvalue']:.2e}" if brain_hit is not None
                    else "  | no brain term in top 50"))
            results.append({
                "list": name, "library": lib,
                "brain_term_rank": int(brain_hit["rank"]) if brain_hit is not None else None,
                "brain_term_adj_p": brain_hit["adj_pvalue"] if brain_hit is not None else np.nan,
                "brain_term_odds_ratio": brain_hit["combined_score"] if brain_hit is not None else np.nan,
            })
            time.sleep(1)
    return pd.DataFrame(results)


def rerun_classifier(pos_symbols, matched_tol_symbols):
    print("\n== Re-running classifier ablation: POS vs tissue-matched TOL ==")
    lab = load_labelled(pos_file=POS_FILE)
    tissue = load_tissue_features()
    joined = lab.join(tissue, on="symbol", how="inner")
    matched_set = set(matched_tol_symbols) | set(pos_symbols)
    matched = joined[joined["symbol"].isin(matched_set)].copy()
    y = matched["label"].values
    print(f"Matched labelled set: n={len(matched)}  POS={int(y.sum())}  TOL={int((y==0).sum())}  "
         f"prevalence={y.mean():.3f}")

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
    sets = {"tissue_only": TISSUE, "level_only": LEVEL, "shape_only": PURE_SHAPE,
           "level+shape": LEVEL + PURE_SHAPE, "level+shape+tissue": LEVEL + PURE_SHAPE + TISSUE}
    rows = []
    print(f"{'feature set':22s} {'ROC-AUC':>16s} {'PR-AUC':>16s}")
    for name, feats in sets.items():
        pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                             scale=True)
        res = cross_validate(pipe, matched[feats], y, cv=cv,
                             scoring=["roc_auc", "average_precision"], n_jobs=-1)
        ra, pa = res["test_roc_auc"].mean(), res["test_average_precision"].mean()
        print(f"{name:22s} {ra:.3f} +/- {res['test_roc_auc'].std():.3f}   "
             f"{pa:.3f} +/- {res['test_average_precision'].std():.3f}")
        rows.append({"feature_set": name, "roc_auc": ra, "pr_auc": pa, "n_feat": len(feats)})
    print(f"(PR-AUC baseline = prevalence = {y.mean():.3f})")
    return pd.DataFrame(rows)


def strict_brain_primary_tol():
    """The strict alternative to nearest-neighbor matching: TOL genes where brain is
    literally the single top-expressed tissue (45 genome-wide, per
    gene_tissue_prominence.csv's definition) -- a much smaller but more literal
    "brain-primary" counterpart to POS than the continuous nearest-neighbor match."""
    t = pd.read_csv("gtex_tissue_specificity.csv")
    t = t.sort_values("top_tissue_tpm", ascending=False).drop_duplicates("symbol", keep="first")
    return set(t.loc[t["top_tissue"].str.startswith("Brain"), "symbol"])


def rerun_classifier_strict(brain_primary_tol_symbols):
    print("\n== Classifier ablation: POS vs STRICT brain-primary TOL (top_tissue==Brain) ==")
    lab = load_labelled(pos_file=POS_FILE)
    tissue = load_tissue_features()
    joined = lab.join(tissue, on="symbol", how="inner")
    pos_symbols = set(joined.loc[joined["label"] == 1, "symbol"])
    keep = pos_symbols | (brain_primary_tol_symbols & set(joined.loc[joined["label"] == 0, "symbol"]))
    strict = joined[joined["symbol"].isin(keep)].copy()
    y = strict["label"].values
    print(f"n={len(strict)}  POS={int(y.sum())}  TOL={int((y==0).sum())}  prevalence={y.mean():.3f}"
         f"  (NOTE: small n, reversed class balance -- treat as exploratory, not a stable estimate)")

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
    sets = {"tissue_only": TISSUE, "level_only": LEVEL, "shape_only": PURE_SHAPE,
           "level+shape": LEVEL + PURE_SHAPE, "level+shape+tissue": LEVEL + PURE_SHAPE + TISSUE}
    rows = []
    print(f"{'feature set':22s} {'ROC-AUC':>16s} {'PR-AUC':>16s}")
    for name, feats in sets.items():
        pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                             scale=True)
        res = cross_validate(pipe, strict[feats], y, cv=cv,
                             scoring=["roc_auc", "average_precision"], n_jobs=-1)
        ra, pa = res["test_roc_auc"].mean(), res["test_average_precision"].mean()
        print(f"{name:22s} {ra:.3f} +/- {res['test_roc_auc'].std():.3f}   "
             f"{pa:.3f} +/- {res['test_average_precision'].std():.3f}")
        rows.append({"feature_set": name, "roc_auc": ra, "pr_auc": pa, "n_feat": len(feats)})
    print(f"(PR-AUC baseline = prevalence = {y.mean():.3f})")
    return pd.DataFrame(rows)


def main():
    pos_t, tol_t, matched_tol = build_tissue_matched_tol()
    enrichr_df = rerun_enrichr(pos_t.index.tolist(), matched_tol.index.tolist())
    enrichr_df.to_csv(SUMMARY_OUT, index=False)
    print(f"\nSaved Enrichr summary -> {SUMMARY_OUT}")

    clf_df = rerun_classifier(pos_t.index.tolist(), matched_tol.index.tolist())
    clf_df.to_csv("tissue_matched_classifier_ablation.csv", index=False)
    print("Saved classifier ablation -> tissue_matched_classifier_ablation.csv")

    strict_tol = strict_brain_primary_tol()
    strict_df = rerun_classifier_strict(strict_tol)
    strict_df.to_csv("tissue_matched_strict_classifier_ablation.csv", index=False)
    print("Saved strict-subset classifier ablation -> tissue_matched_strict_classifier_ablation.csv")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    pivot = enrichr_df.pivot(index="library", columns="list", values="brain_term_adj_p")
    pivot = -np.log10(pivot.clip(lower=1e-300))
    pivot.plot(kind="bar", ax=ax, color=["crimson", "gray"])
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=1, label="p=0.05")
    ax.set_ylabel("-log10(adj p), best brain term")
    ax.set_title("POS vs tissue-matched TOL: brain enrichment")
    ax.legend()
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

    ax = axes[1]
    ax.bar(clf_df["feature_set"], clf_df["roc_auc"], color="steelblue")
    ax.axhline(0.5, ls="--", c="k", lw=1, label="chance")
    ax.set_ylabel("CV ROC-AUC")
    ax.set_title("Classifier ablation: POS vs tissue-matched TOL")
    ax.set_ylim(0.4, 1.0)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.legend()

    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"Saved figure -> {FIG_OUT}")


if __name__ == "__main__":
    main()
