#!/usr/bin/env python
"""Dedicated logistic regression: POS vs the strict brain-primary TOL subset
(TOL genes where brain is their literal single top-expressed tissue -- 45
genome-wide, of which 37 have a valid cerebellum fourparam fit and full
feature set). Follow-up to tissue_matched_reanalysis.py's ablation table --
this fits ONE proper model (not just a CV score per feature-set block) and
reports OOF ROC/confusion + standardized coefficients, so we can see exactly
which features drive the separation in this specific, tissue-matched-by-
construction comparison.

Caveat carried over from tissue_matched_reanalysis.py: n=105 (68 POS + 37 TOL)
with reversed class balance (prevalence 0.648) -- small and unstable, read
directionally, not as a settled estimate.

Run: python logreg_strict_brain_primary.py
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (confusion_matrix, f1_score, precision_recall_curve,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold,
                                     cross_val_predict, cross_validate)

from ablation_dosage_features import LEVEL, PURE_SHAPE
from tissue_matched_reanalysis import strict_brain_primary_tol
from train_dosage_classifier import load_labelled, make_pipeline
from train_tissue_aware_classifier import TISSUE, load_tissue_features

POS_FILE = "data/positive_genes_compiled.txt"
FEATURES = LEVEL + PURE_SHAPE + TISSUE
FIG_OUT = "outputs/figures/logreg_strict_brain_primary_eval.png"
COEF_OUT = "outputs/tables/logreg_strict_brain_primary_coefficients.csv"


def build_data():
    lab = load_labelled(pos_file=POS_FILE)
    tissue = load_tissue_features()
    joined = lab.join(tissue, on="symbol", how="inner")

    brain_primary_genome_wide = strict_brain_primary_tol()  # ALL brain-primary genes genome-wide
    pos_symbols = set(joined.loc[joined["label"] == 1, "symbol"])
    tol_symbols_all = set(lab.loc[lab["label"] == 0, "symbol"])  # full TOL, before tissue join
    brain_primary_tol_all = brain_primary_genome_wide & tol_symbols_all  # restricted to TOL (~45)
    matched_tol = brain_primary_tol_all & set(joined.loc[joined["label"] == 0, "symbol"])

    print(f"Brain-primary TOL genome-wide (45 expected): {len(brain_primary_tol_all)}")
    print(f"...of which have usable cerebellum+tissue features: {len(matched_tol)}")
    missing = brain_primary_tol_all - matched_tol
    if missing:
        print(f"...missing (no valid cerebellum fit): {sorted(missing)}")

    keep = pos_symbols | matched_tol
    df = joined[joined["symbol"].isin(keep)].copy().reset_index(drop=True)
    return df


def main():
    df = build_data()
    X, y = df[FEATURES], df["label"].values
    n_pos, n_tol = int(y.sum()), int((y == 0).sum())
    print(f"\nn={len(y)}  POS={n_pos}  TOL={n_tol}  prevalence={y.mean():.3f}")
    print("NOTE: small n with reversed class balance -- read directionally, not as a "
         "stable estimate.\n")

    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                         scale=True)

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
    res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"], n_jobs=-1)
    print(f"CV ROC-AUC: {res['test_roc_auc'].mean():.3f} +/- {res['test_roc_auc'].std():.3f}")
    print(f"CV PR-AUC : {res['test_average_precision'].mean():.3f} +/- "
         f"{res['test_average_precision'].std():.3f}  (baseline={y.mean():.3f})")

    oof_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    oof = cross_val_predict(pipe, X, y, cv=oof_cv, method="predict_proba", n_jobs=-1)[:, 1]
    roc_auc = roc_auc_score(y, oof)
    prec, rec, thr = precision_recall_curve(y, oof)
    f1s = 2 * prec * rec / (prec + rec + 1e-12)
    best_thr = thr[max(0, np.argmax(f1s[:-1]))]
    yhat = (oof >= best_thr).astype(int)
    cm = confusion_matrix(y, yhat)
    print(f"\nOOF ROC-AUC: {roc_auc:.3f}")
    print(f"OOF F1@thr={best_thr:.3f}: {f1_score(y, yhat):.3f}")
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
             title=f"POS (n={n_pos}) vs brain-primary TOL (n={n_tol})\nOOF ROC")
    ax[0].legend(loc="lower right")

    top = coefs.head(12)[::-1]
    colors = ["crimson" if f in TISSUE else "steelblue" for f in top.index]
    ax[1].barh(top.index, top.values, color=colors)
    ax[1].axvline(0, color="k", lw=0.8)
    ax[1].set(title="Top 12 standardized coefficients (red=tissue feature)",
             xlabel="log-odds per SD")
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"Saved figure -> {FIG_OUT}")


if __name__ == "__main__":
    main()
