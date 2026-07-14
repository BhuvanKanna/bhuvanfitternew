#!/usr/bin/env python
"""Does fitting each gene in its own biologically-relevant tissue (instead of
forcing every gene through cerebellum) unlock a real truncation-index signal?
Baseline to beat: cerebellum-only truncationindex-alone was ROC-AUC 0.528
(~chance) against POS/TOL (see ablation_dosage_features.py).

Adds two scale-free ratios (cv, w_x0) to gene_own_tissue_fourparam.csv (not
computed by generate_own_tissue_fourparam.py itself) and runs the identical
PURE_SHAPE/LEVEL ablation used throughout this repo, on the own-tissue fits
instead of the cerebellum ones.

Run: python compare_own_tissue_vs_cerebellum.py
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate

from ablation_dosage_features import LEVEL, PURE_SHAPE
from train_dosage_classifier import make_pipeline

OWN_TISSUE_TABLE = "gene_own_tissue_fourparam.csv"
OUT = "own_tissue_vs_cerebellum_ablation.csv"

# Cerebellum-baseline numbers already established (ablation_dosage_features.py /
# CLAUDE.md), for direct comparison -- not recomputed here.
CEREBELLUM_BASELINE = {
    "truncationindex_only": 0.528,
    "pure_shape": 0.775,
    "level_only": 0.818,
    "mean_only": 0.824,
}


def main():
    df = pd.read_csv(OWN_TISSUE_TABLE)
    df = df[df["fit_success"] == True].copy()  # noqa: E712
    df["cv"] = df["std"] / df["mean"]
    df["w_x0"] = df["w"] / df["x0"]
    df = df.replace([np.inf, -np.inf], np.nan)

    y_all = df["label"].values
    n_pos, n_tol = int((y_all == 1).sum()), int((y_all == 0).sum())
    print(f"Own-tissue fits with fit_success=True: {len(df)}  (POS={n_pos}, TOL={n_tol})")
    print(f"Tissue breakdown (top 10):")
    print(df["tissue"].value_counts().head(10).to_string())

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
    sets = {
        "truncationindex_only": ["truncationindex"],
        "pure_shape": PURE_SHAPE,
        "level_only": LEVEL,
        "mean_only": ["mean"],
        "level+shape": LEVEL + PURE_SHAPE,
    }
    print(f"\n{'feature set':22s} {'own-tissue ROC-AUC':>20s} {'cerebellum ROC-AUC':>20s}  n_feat")
    rows = []
    for name, feats in sets.items():
        sub = df.dropna(subset=feats)
        y = sub["label"].values
        pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                             scale=True)
        res = cross_validate(pipe, sub[feats], y, cv=cv,
                             scoring=["roc_auc", "average_precision"], n_jobs=-1)
        ra, ras = res["test_roc_auc"].mean(), res["test_roc_auc"].std()
        pa = res["test_average_precision"].mean()
        cereb = CEREBELLUM_BASELINE.get(name)
        cereb_str = f"{cereb:.3f}" if cereb is not None else "n/a"
        print(f"{name:22s} {ra:.3f} +/- {ras:.3f}      {cereb_str:>20s}   {len(feats)}  (n={len(sub)})")
        rows.append({"feature_set": name, "own_tissue_roc_auc": ra, "own_tissue_roc_auc_std": ras,
                     "own_tissue_pr_auc": pa, "cerebellum_roc_auc": cereb, "n": len(sub)})

    out = pd.DataFrame(rows)

    # -- Critical follow-up check: is the raw own-tissue signal genuine gene-level
    # biology, or just a between-tissue artifact (POS/TOL landing in systematically
    # different tissues that happen to have different typical truncationindex/
    # mean/shape values, independent of any real per-gene "ceiling")? Tissue-demean
    # each feature (subtract that tissue's own mean) and re-score -- if the signal
    # collapses, it was between-tissue; if it survives, it's within-tissue (real).
    print("\n== Tissue-demeaned check (does the signal survive WITHIN tissue?) ==")
    demeaned_rows = []
    for name, feats in [("truncationindex_only", ["truncationindex"]),
                        ("mean_only", ["mean"]), ("pure_shape", PURE_SHAPE)]:
        sub = df.dropna(subset=feats).copy()
        for f in feats:
            sub[f + "_dm"] = sub.groupby("tissue")[f].transform(lambda s: s - s.mean())
        dm_feats = [f + "_dm" for f in feats]
        y = sub["label"].values
        pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                             scale=True)
        res = cross_validate(pipe, sub[dm_feats], y, cv=cv, scoring=["roc_auc"], n_jobs=-1)
        raw_auc = out.loc[out.feature_set == name, "own_tissue_roc_auc"].iloc[0]
        dm_auc = res["test_roc_auc"].mean()
        print(f"{name:22s} raw={raw_auc:.3f}  tissue-demeaned={dm_auc:.3f}  "
             f"(delta={dm_auc-raw_auc:+.3f})")
        demeaned_rows.append({"feature_set": name, "raw_roc_auc": raw_auc,
                              "tissue_demeaned_roc_auc": dm_auc})

    demeaned_df = pd.DataFrame(demeaned_rows)
    out = out.merge(demeaned_df[["feature_set", "tissue_demeaned_roc_auc"]],
                    on="feature_set", how="left")
    out.to_csv(OUT, index=False)
    print(f"\nSaved -> {OUT}")

    ti_delta = out.loc[out.feature_set == "truncationindex_only", "own_tissue_roc_auc"].iloc[0] - 0.528
    ti_dm = demeaned_df.loc[demeaned_df.feature_set == "truncationindex_only",
                           "tissue_demeaned_roc_auc"].iloc[0]
    print(f"\ntruncationindex_only: own-tissue(raw) vs cerebellum delta = {ti_delta:+.3f}")
    print(f"truncationindex_only: tissue-demeaned ROC-AUC = {ti_dm:.3f} "
         f"-> {'mostly a between-tissue artifact' if ti_dm < 0.55 else 'signal survives within-tissue'}")


if __name__ == "__main__":
    main()
