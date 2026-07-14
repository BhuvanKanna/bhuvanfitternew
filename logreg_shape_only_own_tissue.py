#!/usr/bin/env python
"""Dedicated shape-only model on the own-tissue fourparam fits: excludes ALL
level/expression-magnitude features (mean, min, max, x0, std, n_obs, A, y0,
maxheight, rightheight) entirely, using only PURE_SHAPE (skew, kurt, cv, w_x0,
truncationindex, ti_fourparam_sigma_dist -- all dimensionless ratios/moments,
none of which encode absolute expression level).

Runs both the raw version and the tissue-demeaned version (subtract each
tissue's own mean from every feature before fitting) -- the raw version can
still be partly inflated by a between-tissue artifact the way truncationindex
alone was (see compare_own_tissue_vs_cerebellum.py), so the demeaned version is
the more trustworthy answer to "does shape alone, with no help from level or
which-tissue-was-used, carry real signal."

Run: python logreg_shape_only_own_tissue.py
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (confusion_matrix, f1_score, precision_recall_curve,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold,
                                     cross_val_predict, cross_validate)

from ablation_dosage_features import PURE_SHAPE
from train_dosage_classifier import make_pipeline

TABLE = "gene_own_tissue_fourparam.csv"
FIG_OUT = "logreg_shape_only_own_tissue_eval.png"
COEF_OUT = "logreg_shape_only_own_tissue_coefficients.csv"


def load_data():
    df = pd.read_csv(TABLE)
    df = df[df["fit_success"] == True].copy()  # noqa: E712
    df["cv"] = df["std"] / df["mean"]
    df["w_x0"] = df["w"] / df["x0"]
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=PURE_SHAPE)
    return df


def run_one(df, feats, label):
    X, y = df[feats], df["label"].values
    n_pos, n_tol = int(y.sum()), int((y == 0).sum())
    print(f"\n== {label} == (n={len(y)}  POS={n_pos}  TOL={n_tol}  prevalence={y.mean():.3f})")

    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                         scale=True)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
    res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"], n_jobs=-1)
    print(f"CV ROC-AUC: {res['test_roc_auc'].mean():.3f} +/- {res['test_roc_auc'].std():.3f}   "
         f"PR-AUC: {res['test_average_precision'].mean():.3f} (baseline={y.mean():.3f})")

    oof_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    oof = cross_val_predict(pipe, X, y, cv=oof_cv, method="predict_proba", n_jobs=-1)[:, 1]
    roc_auc = roc_auc_score(y, oof)
    prec, rec, thr = precision_recall_curve(y, oof)
    f1s = 2 * prec * rec / (prec + rec + 1e-12)
    best_thr = thr[max(0, np.argmax(f1s[:-1]))]
    yhat = (oof >= best_thr).astype(int)
    cm = confusion_matrix(y, yhat)
    print(f"OOF ROC-AUC: {roc_auc:.3f}   F1@thr={f1_score(y, yhat):.3f}")
    print(f"Confusion [[TN FP][FN TP]]:\n{cm}")

    pipe.fit(X, y)
    coefs = pd.Series(pipe.named_steps["clf"].coef_[0], index=feats).sort_values(
        key=np.abs, ascending=False)
    print("Standardized coefficients (log-odds per SD):")
    print(coefs.round(3).to_string())
    return roc_auc, oof, y, coefs


def main():
    df = load_data()

    print("=" * 60)
    print("SHAPE-ONLY MODEL (mean/level fully excluded)")
    print("=" * 60)
    raw_auc, raw_oof, raw_y, raw_coefs = run_one(df, PURE_SHAPE, "RAW (own-tissue, uncontrolled)")

    df_dm = df.copy()
    for f in PURE_SHAPE:
        df_dm[f] = df_dm.groupby("tissue")[f].transform(lambda s: s - s.mean())
    dm_auc, dm_oof, dm_y, dm_coefs = run_one(df_dm, PURE_SHAPE,
                                             "TISSUE-DEMEANED (within-tissue only)")

    raw_coefs.to_csv(COEF_OUT.replace(".csv", "_raw.csv"), header=["coefficient"])
    dm_coefs.to_csv(COEF_OUT.replace(".csv", "_demeaned.csv"), header=["coefficient"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    for oof, y, auc, label, color in [(raw_oof, raw_y, raw_auc, "raw", "crimson"),
                                       (dm_oof, dm_y, dm_auc, "tissue-demeaned", "steelblue")]:
        fpr, tpr, _ = roc_curve(y, oof)
        ax.plot(fpr, tpr, lw=2, color=color, label=f"{label} AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set(xlabel="FPR", ylabel="TPR", title="Shape-only OOF ROC (mean excluded)")
    ax.legend(loc="lower right")

    for ax, coefs, title in [(axes[1], raw_coefs, "Raw: coefficients"),
                             (axes[2], dm_coefs, "Tissue-demeaned: coefficients")]:
        top = coefs[::-1]
        ax.barh(top.index, top.values, color="steelblue")
        ax.axvline(0, color="k", lw=0.8)
        ax.set(title=title, xlabel="log-odds per SD")

    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"\nSaved figure -> {FIG_OUT}")
    print(f"Saved coefficients -> {COEF_OUT.replace('.csv', '_raw.csv')} / _demeaned.csv")

    print(f"\n== Verdict ==")
    print(f"Shape alone (mean/level excluded), raw own-tissue: ROC-AUC={raw_auc:.3f}")
    print(f"Shape alone (mean/level excluded), tissue-demeaned: ROC-AUC={dm_auc:.3f}")
    baseline = df["label"].mean()
    print(f"(prevalence baseline={baseline:.3f}; chance ROC-AUC=0.5 regardless of prevalence)")


if __name__ == "__main__":
    main()
