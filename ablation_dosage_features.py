#!/usr/bin/env python
"""Feature ablation for the dosage classifier: which fourparam columns carry the signal?

Reuses train_dosage_classifier.load_labelled() and compares logistic-regression
CV ROC-AUC / PR-AUC across feature subsets, to answer "which columns are best".
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate

from train_dosage_classifier import load_labelled, make_pipeline

LEVEL = ["mean", "min", "max", "x0", "std", "n_obs", "A", "y0", "maxheight", "rightheight"]
# NOTE: raw `w` is an absolute Gaussian-width parameter -- it scales with expression
# magnitude just like mean/min/max, so it is NOT scale-free. Kept here only to show,
# by contrast, how much it inflates a naive "shape" set (see shape_only_w_contaminated).
SHAPE = ["skew", "kurt", "truncationindex", "ti_fourparam_sigma_dist", "cv", "w_x0", "w"]
# Strictly scale-invariant: unaffected by multiplying all of a gene's expression
# values by a constant. Only dimensionless ratios/moments -- no absolute-unit column.
PURE_SHAPE = ["skew", "kurt", "cv", "w_x0", "truncationindex", "ti_fourparam_sigma_dist"]
SETS = {
    "level_only": LEVEL,
    "shape_only_w_contaminated": SHAPE,
    "pure_shape": PURE_SHAPE,
    "combined": LEVEL + SHAPE,
    "mean_only": ["mean"],
    "mean+truncationindex": ["mean", "truncationindex"],
    "truncationindex_only": ["truncationindex"],
    "top3_level": ["mean", "min", "max"],
}

lab = load_labelled()
y = lab["label"].values
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)

print(f"POS={int(y.sum())}  TOL={int((y==0).sum())}  prevalence={y.mean():.3f}\n")
print(f"{'feature set':22s} {'ROC-AUC':>16s} {'PR-AUC':>16s}  n_feat")
rows = []
for name, feats in SETS.items():
    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=0.01),
                         scale=True)
    res = cross_validate(pipe, lab[feats], y, cv=cv,
                         scoring=["roc_auc", "average_precision"], n_jobs=-1)
    ra, ras = res["test_roc_auc"].mean(), res["test_roc_auc"].std()
    pa, pas = res["test_average_precision"].mean(), res["test_average_precision"].std()
    rows.append((name, ra, pa))
    print(f"{name:22s} {ra:.3f} +/- {ras:.3f}   {pa:.3f} +/- {pas:.3f}   {len(feats)}")

print(f"\n(PR-AUC baseline = prevalence = {y.mean():.3f})")

# C=0.01 above was tuned for the earlier 17-feature model, not for a 6-feature
# pure-shape set -- sweep C dedicated to pure_shape so it isn't handicapped by
# the wrong regularization strength.
print("\n== Dedicated C sweep for pure_shape (does shape alone predict, at its own best C?) ==")
print(f"{'C':>8s} {'ROC-AUC':>16s} {'PR-AUC':>16s}")
best = (None, -1, None)
for C in [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]:
    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=C),
                         scale=True)
    res = cross_validate(pipe, lab[PURE_SHAPE], y, cv=cv,
                         scoring=["roc_auc", "average_precision"], n_jobs=-1)
    ra, ras = res["test_roc_auc"].mean(), res["test_roc_auc"].std()
    pa, pas = res["test_average_precision"].mean(), res["test_average_precision"].std()
    print(f"{C:8.3f} {ra:.3f} +/- {ras:.3f}   {pa:.3f} +/- {pas:.3f}")
    if ra > best[1]:
        best = (C, ra, pa)

print(f"\nBest pure_shape C={best[0]}: ROC-AUC={best[1]:.3f}  PR-AUC={best[2]:.3f}"
      f"  (vs prevalence baseline {y.mean():.3f})")
