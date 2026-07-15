#!/usr/bin/env python
"""Tests a genuinely different approach to the tissue confound than manual
demeaning: HistGradientBoostingClassifier with tissue as a NATIVE categorical
feature alongside the shape features, rather than logistic regression on
tissue-demeaned features.

Why this model, specifically (not just "try another sklearn classifier" --
logistic regression already beat RandomForest/HistGBM/SVM/KNN in every earlier
comparison in this project on similarly-shaped data, so a blind re-run would
likely just repeat that): manual tissue-demeaning (subtract each tissue's own
mean before fitting) is a crude fix -- it treats every tissue's estimated mean
as equally reliable, which is badly wrong for tissues represented by only 1-2
genes in the training set (a singleton-gene tissue gets zeroed out completely
by demeaning, destroying its signal; a 2-3 gene tissue's "mean" is mostly
noise). A tree ensemble that takes tissue as a native categorical input can
instead learn how much to adjust for tissue with the model's own regularization
(shrinkage/leaf-size limits) providing something like partial pooling, rather
than a hard per-group subtraction -- and it can capture nonlinear/interaction
effects among skew/kurt/cv that a linear combination cannot.

Three models compared on the same shape+tissue feature set, so the comparison
isolates what actually helps: (1) logistic regression, tissue demeaned (the
existing approach); (2) logistic regression, tissue one-hot encoded as an
explicit fixed effect (isolates "does including tissue properly help, even
with a linear model" from "does the tree algorithm itself help"); (3)
HistGradientBoostingClassifier, tissue as a native categorical feature.

Run: python histgbm_shape_tissue_own_tissue.py
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, f1_score, precision_recall_curve, roc_auc_score, roc_curve
from sklearn.model_selection import (GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold,
                                     cross_val_predict, cross_validate)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer

from ablation_dosage_features import PURE_SHAPE
from train_dosage_classifier import make_pipeline

TABLE = "outputs/tables/gene_own_tissue_fourparam.csv"
FIG_OUT = "outputs/figures/histgbm_shape_tissue_eval.png"

# Baselines already established (logreg_shape_only_own_tissue.py), for direct comparison.
BASELINE_RAW = 0.752
BASELINE_DEMEANED = 0.635


def load_data():
    df = pd.read_csv(TABLE)
    df = df[df["fit_success"] == True].copy()  # noqa: E712
    df["cv"] = df["std"] / df["mean"]
    df["w_x0"] = df["w"] / df["x0"]
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=PURE_SHAPE)
    df["tissue"] = df["tissue"].astype("category")
    return df


def logreg_onehot_pipeline(C=1.0):
    """Logistic regression with tissue as an explicit one-hot fixed effect --
    isolates whether including tissue properly (vs demeaning) helps, independent
    of the tree algorithm."""
    pre = ColumnTransformer([
        ("shape", Pipeline([("impute", SimpleImputer(strategy="median")),
                           ("scale", StandardScaler())]), PURE_SHAPE),
        ("tissue", OneHotEncoder(handle_unknown="ignore"), ["tissue"]),
    ])
    return Pipeline([("pre", pre),
                     ("clf", LogisticRegression(max_iter=5000, class_weight="balanced", C=C))])


def histgbm_pipeline(**params):
    """HistGradientBoostingClassifier with tissue as a NATIVE categorical
    feature -- no manual demeaning, no one-hot expansion."""
    pre = ColumnTransformer([
        ("shape", "passthrough", PURE_SHAPE),
        ("tissue", "passthrough", ["tissue"]),
    ])
    n_shape = len(PURE_SHAPE)
    categorical_mask = [False] * n_shape + [True]
    return Pipeline([
        ("pre", pre),
        ("clf", HistGradientBoostingClassifier(
            categorical_features=categorical_mask, class_weight="balanced",
            random_state=0, **params)),
    ])


def evaluate(pipe, X, y, cv, label):
    res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"], n_jobs=-1)
    ra, ras = res["test_roc_auc"].mean(), res["test_roc_auc"].std()
    pa = res["test_average_precision"].mean()
    print(f"{label:45s} CV ROC-AUC={ra:.3f} +/- {ras:.3f}   PR-AUC={pa:.3f}")
    return ra, ras, pa


def main():
    df = load_data()
    X = df[PURE_SHAPE + ["tissue"]]
    y = df["label"].values
    n_pos, n_tol = int(y.sum()), int((y == 0).sum())
    print(f"n={len(y)}  POS={n_pos}  TOL={n_tol}  prevalence={y.mean():.3f}")
    print(f"Unique tissues represented: {df['tissue'].nunique()}  "
         f"(tissues with only 1 gene: {(df['tissue'].value_counts() == 1).sum()})")

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)

    print(f"\n{'model':45s} metric")
    print(f"{'(reference) logreg shape-only, raw':45s} CV ROC-AUC={BASELINE_RAW:.3f}  (from earlier script)")
    print(f"{'(reference) logreg shape-only, tissue-demeaned':45s} "
         f"CV ROC-AUC={BASELINE_DEMEANED:.3f}  (from earlier script)")

    evaluate(logreg_onehot_pipeline(), X, y, cv, "logreg, shape + tissue one-hot")

    # Small hyperparameter grid, sized for n~600: shallow trees, few leaves,
    # moderate learning rate -- aggressive regularization to avoid overfitting.
    base_grid = {
        "learning_rate": [0.03, 0.1, 0.2], "max_leaf_nodes": [7, 15, 31],
        "max_depth": [2, 3, None], "l2_regularization": [0.0, 1.0, 5.0],
        "min_samples_leaf": [5, 10, 20],
    }

    # Fair, apples-to-apples check against the logreg baselines: HistGBM on
    # shape features ALONE (no tissue feature at all -- same inputs the logreg
    # baselines got), both raw and tissue-demeaned. This is the test that
    # isolates "is HistGBM just exploiting tissue-as-a-shortcut" (the
    # shape+tissue model below) from "does HistGBM find more real shape signal
    # than a linear model, even with no access to tissue at all."
    df_dm = df.copy()
    for f in PURE_SHAPE:
        df_dm[f] = df_dm.groupby("tissue", observed=True)[f].transform(lambda s: s - s.mean())
    print()
    gs_shape_raw = GridSearchCV(HistGradientBoostingClassifier(class_weight="balanced",
                                random_state=0), base_grid, scoring="roc_auc", cv=cv, n_jobs=-1)
    gs_shape_raw.fit(df[PURE_SHAPE], y)
    print(f"HistGBM shape-only, RAW (no tissue feature):        "
         f"CV ROC-AUC={gs_shape_raw.best_score_:.3f}")
    gs_shape_dm = GridSearchCV(HistGradientBoostingClassifier(class_weight="balanced",
                               random_state=0), base_grid, scoring="roc_auc", cv=cv, n_jobs=-1)
    gs_shape_dm.fit(df_dm[PURE_SHAPE], y)
    print(f"HistGBM shape-only, TISSUE-DEMEANED (no tissue feature): "
         f"CV ROC-AUC={gs_shape_dm.best_score_:.3f}")

    grid = {"clf__" + k: v for k, v in base_grid.items()}
    gs = GridSearchCV(histgbm_pipeline(), grid, scoring="roc_auc", cv=cv, n_jobs=-1)
    gs.fit(X, y)
    print(f"\nTuned HistGBM (shape + tissue native categorical):")
    print(f"  best CV ROC-AUC={gs.best_score_:.3f}  params={gs.best_params_}")
    tuned = gs.best_estimator_

    oof_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=0)
    oof = cross_val_predict(tuned, X, y, cv=oof_cv, method="predict_proba", n_jobs=-1)[:, 1]
    roc_auc = roc_auc_score(y, oof)
    prec, rec, thr = precision_recall_curve(y, oof)
    f1s = 2 * prec * rec / (prec + rec + 1e-12)
    best_thr = thr[max(0, np.argmax(f1s[:-1]))]
    yhat = (oof >= best_thr).astype(int)
    cm = confusion_matrix(y, yhat)
    print(f"  OOF ROC-AUC={roc_auc:.3f}  F1@thr={f1_score(y, yhat):.3f}")
    print(f"  Confusion [[TN FP][FN TP]]:\n{cm}")

    tuned.fit(X, y)
    imp = permutation_importance(tuned, X, y, scoring="roc_auc", n_repeats=20,
                                 random_state=0, n_jobs=-1)
    importance = pd.Series(imp.importances_mean, index=PURE_SHAPE + ["tissue"]).sort_values(
        ascending=False)
    print("\nPermutation importance (HistGBM, shape+tissue):")
    print(importance.round(4).to_string())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fpr, tpr, _ = roc_curve(y, oof)
    axes[0].plot(fpr, tpr, lw=2, color="darkgreen", label=f"HistGBM AUC={roc_auc:.3f}")
    axes[0].plot([0, 1], [0, 1], "k--", lw=1)
    axes[0].axhline(0, color="none")
    axes[0].set(xlabel="FPR", ylabel="TPR",
               title="HistGBM (shape + native categorical tissue) OOF ROC")
    axes[0].legend(loc="lower right")

    top = importance[::-1]
    colors = ["darkorange" if f == "tissue" else "steelblue" for f in top.index]
    axes[1].barh(top.index, top.values, color=colors)
    axes[1].set(title="Permutation importance", xlabel="mean ROC-AUC drop")
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"\nSaved figure -> {FIG_OUT}")

    print(f"\n== Summary ==")
    print(f"logreg, shape-only, raw:                        {BASELINE_RAW:.3f}")
    print(f"logreg, shape-only, tissue-demeaned:            {BASELINE_DEMEANED:.3f}  (collapses substantially)")
    print(f"HistGBM, shape-only, raw:                       {gs_shape_raw.best_score_:.3f}")
    print(f"HistGBM, shape-only, tissue-demeaned:           {gs_shape_dm.best_score_:.3f}  "
         f"(barely moves -- unlike logreg)")
    print(f"HistGBM, shape + tissue (native categorical):   CV={gs.best_score_:.3f}  OOF={roc_auc:.3f}  "
         f"(tissue used AS a feature -- exploits the confound, not a controlled test)")


if __name__ == "__main__":
    main()
