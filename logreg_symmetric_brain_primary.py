#!/usr/bin/env python
"""Symmetric fix to logreg_strict_brain_primary.py's construction artifact: that
script restricted TOL to brain-primary genes but left POS unrestricted, so
brain_is_top was true by definition for 100% of TOL but only ~28% of POS --
the model was substantially re-deriving its own selection rule rather than
finding real signal.

This version restricts BOTH sides to brain-primary genes (POS: 19 of 70 with
usable data; TOL: 37 of 45 genome-wide with usable data) so brain_is_top is
true for 100% of both groups -- constant within the dataset, contributing
nothing to discrimination, exactly the intended fix. Whatever signal remains
must come from features OTHER than "is brain my top tissue".

Caveat: n=56 (19 POS + 37 TOL) is very small -- read directionally only.

Run: python logreg_symmetric_brain_primary.py
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (confusion_matrix, f1_score, precision_recall_curve,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import (LeaveOneOut, RepeatedStratifiedKFold,
                                     cross_val_predict, cross_validate)

from ablation_dosage_features import LEVEL, PURE_SHAPE
from train_dosage_classifier import load_labelled, make_pipeline
from train_tissue_aware_classifier import TISSUE, load_tissue_features

POS_FILE = "positive_genes_compiled.txt"
# brain_is_top dropped: constant (=1) for every gene in this dataset by construction,
# zero variance, cannot discriminate -- keeping it would just add noise/collinearity.
TISSUE_NO_TOP = [f for f in TISSUE if f != "brain_is_top"]
FEATURES = LEVEL + PURE_SHAPE + TISSUE_NO_TOP
FIG_OUT = "logreg_symmetric_brain_primary_eval.png"
COEF_OUT = "logreg_symmetric_brain_primary_coefficients.csv"


def build_data():
    lab = load_labelled(pos_file=POS_FILE)
    tissue = load_tissue_features()
    joined = lab.join(tissue, on="symbol", how="inner")
    df = joined[joined["brain_is_top"] == 1].copy().reset_index(drop=True)
    return df


def main():
    df = build_data()
    X, y = df[FEATURES], df["label"].values
    n_pos, n_tol = int(y.sum()), int((y == 0).sum())
    print(f"Symmetric brain-primary set: n={len(y)}  POS={n_pos}  TOL={n_tol}  "
         f"prevalence={y.mean():.3f}")
    print("brain_is_top dropped (constant=1 for all rows by construction).")
    print("NOTE: very small n -- read directionally only, not as a stable estimate.\n")

    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                         scale=True)

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=20, random_state=0)
    res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"], n_jobs=-1)
    print(f"CV ROC-AUC: {res['test_roc_auc'].mean():.3f} +/- {res['test_roc_auc'].std():.3f}")
    print(f"CV PR-AUC : {res['test_average_precision'].mean():.3f} +/- "
         f"{res['test_average_precision'].std():.3f}  (baseline={y.mean():.3f})")

    # n=56 is too small for a stable 5-fold OOF split -- leave-one-out uses every
    # gene as its own held-out test case, the most data-efficient option here.
    loo = LeaveOneOut()
    oof = cross_val_predict(pipe, X, y, cv=loo, method="predict_proba", n_jobs=-1)[:, 1]
    roc_auc = roc_auc_score(y, oof)
    prec, rec, thr = precision_recall_curve(y, oof)
    f1s = 2 * prec * rec / (prec + rec + 1e-12)
    best_thr = thr[max(0, np.argmax(f1s[:-1]))]
    yhat = (oof >= best_thr).astype(int)
    cm = confusion_matrix(y, yhat)
    print(f"\nLOO OOF ROC-AUC: {roc_auc:.3f}")
    print(f"LOO OOF F1@thr={best_thr:.3f}: {f1_score(y, yhat):.3f}")
    print(f"Confusion [[TN FP][FN TP]]:\n{cm}")

    pipe.fit(X, y)
    coefs = pd.Series(pipe.named_steps["clf"].coef_[0], index=FEATURES)
    coefs = coefs.sort_values(key=np.abs, ascending=False)
    print("\n== Standardized logistic coefficients (log-odds per SD) ==")
    print(coefs.round(3).to_string())
    coefs.to_csv(COEF_OUT, header=["coefficient"])
    print(f"\nSaved coefficients -> {COEF_OUT}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    fpr, tpr, _ = roc_curve(y, oof)
    ax[0].plot(fpr, tpr, lw=2, color="purple", label=f"AUC={roc_auc:.3f}")
    ax[0].plot([0, 1], [0, 1], "k--", lw=1)
    ax[0].set(xlabel="FPR", ylabel="TPR",
             title=f"POS (n={n_pos}) vs TOL (n={n_tol}), both brain-primary\nLOO OOF ROC")
    ax[0].legend(loc="lower right")

    top = coefs.head(12)[::-1]
    colors = ["crimson" if f in TISSUE_NO_TOP else "steelblue" for f in top.index]
    ax[1].barh(top.index, top.values, color=colors)
    ax[1].axvline(0, color="k", lw=0.8)
    ax[1].set(title="Top 12 standardized coefficients (red=tissue feature)",
             xlabel="log-odds per SD")
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"Saved figure -> {FIG_OUT}")


if __name__ == "__main__":
    main()
