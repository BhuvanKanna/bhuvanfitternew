#!/usr/bin/env python
"""Train a 0-1 dosage-sensitivity classifier from the cerebellum fourparam table.

Positive set  = ``positive_genes.txt``            (curated OE / dosage-sensitive, "POS")
Negative set  = ``positiveANDnegativeControlGenes.csv`` (duplication-tolerant, "TOL")

Each gene is described only by its per-gene expression-distribution features (the
columns of ``cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv`` plus two
scale-free ratios). We map the human symbols in POS/TOL onto the table's ENSG rows
via the GTEx ``Name``/``Description`` id columns (same bridge section 5-7 of the
cerebellum notebook uses), keep the canonical fit_success row per symbol, then:

  1. compare Logistic / RandomForest / HistGBM by repeated-stratified-CV ROC-AUC & PR-AUC,
  2. tune the winning family,
  3. report out-of-fold ROC/PR/confusion + feature importance,
  4. refit on all labelled genes and score EVERY fit_success gene genome-wide,
  5. save the fitted pipeline (joblib) and a genome-wide predictions CSV.

IMPORTANT on "accuracy": the classes are ~11:1 (42 POS vs 472 TOL), so a model that
always predicts TOL already scores ~92% accuracy while being useless. ROC-AUC and
PR-AUC (average precision) are the honest headline metrics; accuracy / balanced
accuracy are reported at an F1-optimal threshold for completeness only.

Run:  python train_dosage_classifier.py [--no-genome] [--seed N]
"""
from __future__ import annotations

import argparse
import json
import warnings

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedStratifiedKFold,
    cross_val_predict,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=UserWarning)

TABLE = "cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv"
SOURCE = "cerebellumlog2.csv"
POS_FILE = "positive_genes.txt"
TOL_FILE = "positiveANDnegativeControlGenes.csv"

MODEL_OUT = "dosage_classifier.joblib"
PRED_OUT = "cerebellum_dosage_classifier_pred.csv"
METRICS_OUT = "dosage_classifier_metrics.json"
FIG_OUT = "dosage_classifier_eval.png"

# Table-only features (no re-histogramming). Scale-free ratios cv & w_x0 are derived.
BASE_FEATURES = [
    "y0", "A", "x0", "w", "ti_fourparam_sigma_dist", "truncationindex",
    "min", "max", "mean", "std", "skew", "kurt", "maxheight", "rightheight", "n_obs",
]
DERIVED = ["cv", "w_x0"]
FEATURES = BASE_FEATURES + DERIVED


def symbol_map():
    """Version-stripped ENSG -> gene symbol, from the GTEx source's two id columns."""
    idm = pd.read_csv(SOURCE, usecols=["Name", "Description"])
    return dict(zip(idm["Name"].str.split(".").str[0], idm["Description"]))


def load_labelled(pos_file=POS_FILE):
    """Return (X, y, meta_df) for POS/TOL genes with one canonical fit_success row each.

    `pos_file` defaults to the 44-gene curated OMIM/G2P list; pass
    "positive_genes_compiled.txt" for the 70-gene list expanded with
    grant_genes.csv's mcOE-sensitive genes (see train_tissue_aware_classifier.py)."""
    sym = symbol_map()
    df = pd.read_csv(TABLE)
    df = df[df["fit_success"] == True].copy()  # noqa: E712
    df["symbol"] = df["gene"].str.split(".").str[0].map(sym)
    df = add_derived(df)

    pos = {l.strip() for l in open(pos_file) if l.strip()}
    tol_raw = pd.read_csv(TOL_FILE, header=None)[0].dropna().astype(str).str.strip()
    tol = set(tol_raw) - {"geneDUPtol"}
    tol -= pos  # POS wins any (there is currently no overlap)

    df["label"] = np.where(df["symbol"].isin(pos), 1,
                           np.where(df["symbol"].isin(tol), 0, -1))
    lab = df[df["label"] >= 0].copy()
    # One row per symbol: the fit_success ENSG with the most observations.
    lab = (lab.sort_values("n_obs", ascending=False)
              .drop_duplicates("symbol", keep="first")
              .reset_index(drop=True))
    return lab


def add_derived(df):
    out = df.copy()
    out["cv"] = out["std"] / out["mean"]
    out["w_x0"] = out["w"] / out["x0"]
    return out.replace([np.inf, -np.inf], np.nan)


def make_pipeline(estimator, scale):
    steps = [("impute", SimpleImputer(strategy="median"))]
    if scale:
        steps.append(("scale", StandardScaler()))
    steps.append(("clf", estimator))
    return Pipeline(steps)


def candidate_models(seed):
    return {
        "logreg": (make_pipeline(
            LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
            scale=True), True),
        "rf": (make_pipeline(
            RandomForestClassifier(n_estimators=400, class_weight="balanced_subsample",
                                   random_state=seed, n_jobs=-1),
            scale=False), False),
        "histgbm": (make_pipeline(
            HistGradientBoostingClassifier(random_state=seed,
                                           class_weight="balanced"),
            scale=False), False),
    }


def cv_scores(pipe, X, y, cv):
    """Mean +/- std ROC-AUC and PR-AUC across all CV folds."""
    from sklearn.model_selection import cross_validate
    res = cross_validate(pipe, X, y, cv=cv,
                         scoring=["roc_auc", "average_precision"], n_jobs=-1)
    return (res["test_roc_auc"].mean(), res["test_roc_auc"].std(),
            res["test_average_precision"].mean(), res["test_average_precision"].std())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-genome", action="store_true",
                    help="skip the genome-wide scoring / CSV emit")
    args = ap.parse_args()

    lab = load_labelled()
    X = lab[FEATURES]
    y = lab["label"].values
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    prevalence = n_pos / len(y)
    print(f"Labelled genes: {len(y)}  (POS={n_pos}, TOL={n_neg}, "
          f"prevalence={prevalence:.3f}, all-negative accuracy={1 - prevalence:.3f})")
    print(f"Features ({len(FEATURES)}): {FEATURES}")
    print(f"NaNs per feature: {X.isna().sum().to_dict()}\n")

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=args.seed)

    print("== Model comparison (repeated 5-fold x10 CV) ==")
    print(f"{'model':10s} {'ROC-AUC':>16s} {'PR-AUC':>16s}")
    scored = {}
    for name, (pipe, _scale) in candidate_models(args.seed).items():
        ra, ras, pa, pas = cv_scores(pipe, X, y, cv)
        scored[name] = ra
        print(f"{name:10s} {ra:.3f} +/- {ras:.3f}   {pa:.3f} +/- {pas:.3f}")
    print(f"(PR-AUC baseline = prevalence = {prevalence:.3f})\n")

    best = max(scored, key=scored.get)
    print(f"Best family by ROC-AUC: {best}\n")

    # ---- Hyperparameter tuning of the winning family -------------------------
    base_pipe, _scale = candidate_models(args.seed)[best]
    grids = {
        "logreg": {"clf__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]},
        "rf": {"clf__n_estimators": [300, 600],
               "clf__max_depth": [None, 4, 8],
               "clf__min_samples_leaf": [1, 3, 5]},
        "histgbm": {"clf__learning_rate": [0.03, 0.1, 0.2],
                    "clf__max_depth": [None, 3, 5],
                    "clf__max_leaf_nodes": [15, 31],
                    "clf__l2_regularization": [0.0, 1.0]},
    }
    gs = GridSearchCV(base_pipe, grids[best], scoring="roc_auc", cv=cv, n_jobs=-1)
    gs.fit(X, y)
    print(f"Tuned {best}: best CV ROC-AUC={gs.best_score_:.3f}  params={gs.best_params_}\n")
    tuned = gs.best_estimator_

    # ---- Out-of-fold predictions for honest metrics --------------------------
    oof_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=args.seed)
    oof = cross_val_predict(tuned, X, y, cv=oof_cv, method="predict_proba", n_jobs=-1)[:, 1]
    roc_auc = roc_auc_score(y, oof)
    pr_auc = average_precision_score(y, oof)

    # F1-optimal threshold from the OOF scores.
    prec, rec, thr = precision_recall_curve(y, oof)
    f1s = 2 * prec * rec / (prec + rec + 1e-12)
    best_thr = thr[max(0, np.argmax(f1s[:-1]))]
    yhat = (oof >= best_thr).astype(int)
    cm = confusion_matrix(y, yhat)
    bal_acc = balanced_accuracy_score(y, yhat)
    acc = (yhat == y).mean()
    f1 = f1_score(y, yhat)

    print("== Out-of-fold performance (single 5-fold split) ==")
    print(f"ROC-AUC          : {roc_auc:.3f}")
    print(f"PR-AUC (avg prec): {pr_auc:.3f}   (baseline {prevalence:.3f})")
    print(f"F1-optimal thr   : {best_thr:.3f}")
    print(f"Accuracy@thr     : {acc:.3f}   (all-negative baseline {1 - prevalence:.3f})")
    print(f"Balanced acc@thr : {bal_acc:.3f}")
    print(f"F1@thr           : {f1:.3f}")
    print(f"Confusion [ [TN FP] [FN TP] ]:\n{cm}\n")

    # ---- Feature importance --------------------------------------------------
    tuned.fit(X, y)
    imp = permutation_importance(tuned, X, y, scoring="roc_auc",
                                 n_repeats=20, random_state=args.seed, n_jobs=-1)
    importance = (pd.Series(imp.importances_mean, index=FEATURES)
                  .sort_values(ascending=False))
    print("== Permutation importance (drop in ROC-AUC when shuffled) ==")
    print(importance.round(4).to_string())
    if best == "logreg":
        coef = tuned.named_steps["clf"].coef_[0]
        print("\nStandardized logistic coefficients (log-odds per SD):")
        print(pd.Series(coef, index=FEATURES).sort_values(key=np.abs, ascending=False)
              .round(3).to_string())
    print()

    # ---- Diagnostic figure ---------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    fpr, tpr, _ = roc_curve(y, oof)
    ax[0].plot(fpr, tpr, lw=2, label=f"AUC={roc_auc:.3f}")
    ax[0].plot([0, 1], [0, 1], "k--", lw=1)
    ax[0].set(xlabel="FPR", ylabel="TPR", title="OOF ROC")
    ax[0].legend(loc="lower right")
    ax[1].plot(rec, prec, lw=2, label=f"AP={pr_auc:.3f}")
    ax[1].axhline(prevalence, ls="--", c="k", lw=1, label=f"baseline={prevalence:.3f}")
    ax[1].set(xlabel="Recall", ylabel="Precision", title="OOF Precision-Recall")
    ax[1].legend(loc="upper right")
    top = importance.head(12)[::-1]
    ax[2].barh(top.index, top.values, color="steelblue")
    ax[2].set(title="Permutation importance (top 12)", xlabel="mean ROC-AUC drop")
    fig.suptitle(f"Dosage classifier ({best}) — POS={n_pos} vs TOL={n_neg}", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"Saved figure -> {FIG_OUT}")

    # ---- Persist model + metrics --------------------------------------------
    import joblib
    joblib.dump({"pipeline": tuned, "features": FEATURES, "best_family": best,
                 "f1_threshold": float(best_thr)}, MODEL_OUT)
    metrics = {
        "n_pos": n_pos, "n_neg": n_neg, "prevalence": prevalence,
        "model": best, "best_params": gs.best_params_,
        "cv_roc_auc_tuned": float(gs.best_score_),
        "oof_roc_auc": float(roc_auc), "oof_pr_auc": float(pr_auc),
        "f1_threshold": float(best_thr), "oof_accuracy": float(acc),
        "oof_balanced_accuracy": float(bal_acc), "oof_f1": float(f1),
        "confusion_matrix": cm.tolist(),
        "model_comparison_roc_auc": {k: float(v) for k, v in scored.items()},
        "permutation_importance": importance.round(5).to_dict(),
    }
    json.dump(metrics, open(METRICS_OUT, "w"), indent=2)
    print(f"Saved model  -> {MODEL_OUT}")
    print(f"Saved metrics-> {METRICS_OUT}")

    # ---- Genome-wide scoring -------------------------------------------------
    if not args.no_genome:
        sym = symbol_map()
        allg = pd.read_csv(TABLE)
        allg = allg[allg["fit_success"] == True].copy()  # noqa: E712
        allg = add_derived(allg)
        allg["symbol"] = allg["gene"].str.split(".").str[0].map(sym)
        score = tuned.predict_proba(allg[FEATURES])[:, 1]
        allg["dosage_score"] = score
        posset = {l.strip() for l in open(POS_FILE) if l.strip()}
        tol_raw = pd.read_csv(TOL_FILE, header=None)[0].dropna().astype(str).str.strip()
        tolset = set(tol_raw) - {"geneDUPtol"}
        allg["label"] = np.where(allg["symbol"].isin(posset), "POS",
                                 np.where(allg["symbol"].isin(tolset), "TOL", ""))
        out = allg[["gene", "symbol", "dosage_score", "label", "n_obs"]].copy()
        out = out.sort_values("dosage_score", ascending=False).reset_index(drop=True)
        out.to_csv(PRED_OUT, index=False)
        print(f"\nGenome-wide scores -> {PRED_OUT}  ({len(out)} fit_success genes)")
        print("Top 15 predicted dosage-sensitive genes:")
        print(out.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
