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
SHAPE = ["skew", "kurt", "truncationindex", "ti_fourparam_sigma_dist", "cv", "w_x0", "w"]
SETS = {
    "level_only": LEVEL,
    "shape_only": SHAPE,
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
