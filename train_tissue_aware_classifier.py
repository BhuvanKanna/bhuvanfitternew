#!/usr/bin/env python
"""Tissue-aware dosage-sensitivity classifier: tests the professor's hypothesis
that OE-tolerant genes are less prominently expressed in brain (relative to
other tissues) than OE-sensitive genes, and folds the result into a final 0-1
dosage-sensitivity model.

Joins the existing cerebellum-fourparam LEVEL/PURE_SHAPE features (now using the
70-gene expanded positive list, positive_genes_compiled.txt = positive_genes.txt
+ grant_genes.csv's 27 mcOE-sensitive genes) with genome-wide GTEx multi-tissue
features (gtex_tissue_specificity.csv) by symbol, then:
  1. ablation: tissue_only / level_only / shape_only / level+shape / level+shape+tissue
     -- the direct test of whether brain-tissue-specificity adds signal beyond
     the existing level/shape features, or is redundant with them.
  2. model comparison + tuning on the winning feature combination.
  3. OOF ROC/PR/confusion, permutation importance, diagnostic figure.
  4. genome-wide 0-1 scoring on every fit_success cerebellum gene with a tissue match.

Run: python train_tissue_aware_classifier.py
"""
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, confusion_matrix, f1_score,
                             precision_recall_curve, roc_auc_score, roc_curve)
from sklearn.model_selection import (GridSearchCV, RepeatedStratifiedKFold,
                                     cross_val_predict, cross_validate)

warnings.filterwarnings("ignore", category=UserWarning)

from ablation_dosage_features import LEVEL, PURE_SHAPE
from train_dosage_classifier import (TABLE, add_derived, load_labelled, make_pipeline,
                                     symbol_map)

POS_FILE = "data/positive_genes_compiled.txt"
TISSUE_TABLE = "outputs/tables/gtex_tissue_specificity.csv"

MODEL_OUT = "outputs/models/dosage_classifier_tissue_aware.joblib"
PRED_OUT = "outputs/tables/dosage_score_tissue_aware.csv"
METRICS_OUT = "outputs/models/dosage_classifier_tissue_aware_metrics.json"
FIG_OUT = "outputs/figures/dosage_classifier_tissue_aware_eval.png"

TISSUE = ["brain_cerebellum", "brain_mean", "brain_max", "nonbrain_mean",
         "nonbrain_max", "log2fc_brain_vs_nonbrain", "tau", "brain_is_top"]


def load_tissue_features():
    t = pd.read_csv(TISSUE_TABLE)
    t["brain_is_top"] = t["top_tissue"].str.startswith("Brain").astype(int)
    # one row per symbol: highest top_tissue_tpm (arbitrary canonical tie-break)
    t = t.sort_values("top_tissue_tpm", ascending=False).drop_duplicates("symbol", keep="first")
    return t.set_index("symbol")[TISSUE]


def build_joined_table():
    lab = load_labelled(pos_file=POS_FILE)
    tissue = load_tissue_features()
    joined = lab.join(tissue, on="symbol", how="inner")
    return joined


def cv_scores(X, y, cv, C=1.0):
    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=C),
                         scale=True)
    res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"], n_jobs=-1)
    return res["test_roc_auc"].mean(), res["test_roc_auc"].std(), \
           res["test_average_precision"].mean(), res["test_average_precision"].std()


def run_ablation(df, y, cv):
    sets = {
        "tissue_only": TISSUE,
        "level_only": LEVEL,
        "shape_only": PURE_SHAPE,
        "level+shape": LEVEL + PURE_SHAPE,
        "level+shape+tissue": LEVEL + PURE_SHAPE + TISSUE,
    }
    print(f"\n{'feature set':22s} {'ROC-AUC':>16s} {'PR-AUC':>16s}  n_feat")
    results = {}
    for name, feats in sets.items():
        ra, ras, pa, pas = cv_scores(df[feats], y, cv)
        results[name] = ra
        print(f"{name:22s} {ra:.3f} +/- {ras:.3f}   {pa:.3f} +/- {pas:.3f}   {len(feats)}")
    print(f"(PR-AUC baseline = prevalence = {y.mean():.3f})")
    return results


def main():
    df = build_joined_table()
    y = df["label"].values
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    print(f"Joined labelled genes: {len(y)} (POS={n_pos}, TOL={n_neg}, "
         f"prevalence={y.mean():.3f})")
    print(f"Lost to missing tissue match: "
         f"{len(load_labelled(pos_file=POS_FILE)) - len(df)} genes")

    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)

    print("\n== Ablation: does brain-tissue-specificity add signal? ==")
    ablation = run_ablation(df, y, cv)

    best_combo_name = max(ablation, key=ablation.get)
    FEATURES = {
        "tissue_only": TISSUE, "level_only": LEVEL, "shape_only": PURE_SHAPE,
        "level+shape": LEVEL + PURE_SHAPE, "level+shape+tissue": LEVEL + PURE_SHAPE + TISSUE,
    }[best_combo_name]
    print(f"\nBest feature combination: {best_combo_name} ({len(FEATURES)} features)")

    # ---- Model comparison on the winning feature combination ------------------
    X = df[FEATURES]
    print("\n== Model comparison (repeated 5-fold x10 CV) ==")
    candidates = {
        "logreg": make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced"),
                                scale=True),
        "rf": make_pipeline(RandomForestClassifier(n_estimators=400,
                            class_weight="balanced_subsample", random_state=0, n_jobs=-1),
                            scale=False),
        "histgbm": make_pipeline(HistGradientBoostingClassifier(random_state=0,
                                 class_weight="balanced"), scale=False),
    }
    scored = {}
    for name, pipe in candidates.items():
        res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"],
                             n_jobs=-1)
        ra = res["test_roc_auc"].mean()
        scored[name] = ra
        print(f"{name:10s} ROC-AUC={ra:.3f} +/- {res['test_roc_auc'].std():.3f}  "
             f"PR-AUC={res['test_average_precision'].mean():.3f}")
    best = max(scored, key=scored.get)
    print(f"Best model family: {best}")

    # ---- Tune the winner -------------------------------------------------------
    grids = {
        "logreg": {"clf__C": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]},
        "rf": {"clf__n_estimators": [300, 600], "clf__max_depth": [None, 4, 8],
              "clf__min_samples_leaf": [1, 3, 5]},
        "histgbm": {"clf__learning_rate": [0.03, 0.1, 0.2], "clf__max_depth": [None, 3, 5],
                   "clf__max_leaf_nodes": [15, 31], "clf__l2_regularization": [0.0, 1.0]},
    }
    gs = GridSearchCV(candidates[best], grids[best], scoring="roc_auc", cv=cv, n_jobs=-1)
    gs.fit(X, y)
    print(f"\nTuned {best}: best CV ROC-AUC={gs.best_score_:.3f}  params={gs.best_params_}")
    tuned = gs.best_estimator_

    # ---- OOF metrics ------------------------------------------------------------
    oof_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=0)
    oof = cross_val_predict(tuned, X, y, cv=oof_cv, method="predict_proba", n_jobs=-1)[:, 1]
    roc_auc = roc_auc_score(y, oof)
    pr_auc = average_precision_score(y, oof)
    prec, rec, thr = precision_recall_curve(y, oof)
    f1s = 2 * prec * rec / (prec + rec + 1e-12)
    best_thr = thr[max(0, np.argmax(f1s[:-1]))]
    yhat = (oof >= best_thr).astype(int)
    cm = confusion_matrix(y, yhat)
    print(f"\n== Out-of-fold performance ==\nROC-AUC={roc_auc:.3f}  PR-AUC={pr_auc:.3f}  "
         f"F1@thr={f1_score(y, yhat):.3f}\nConfusion [[TN FP][FN TP]]:\n{cm}")

    # ---- Feature importance ------------------------------------------------------
    tuned.fit(X, y)
    imp = permutation_importance(tuned, X, y, scoring="roc_auc", n_repeats=20,
                                 random_state=0, n_jobs=-1)
    importance = pd.Series(imp.importances_mean, index=FEATURES).sort_values(ascending=False)
    print("\n== Permutation importance ==")
    print(importance.round(4).to_string())

    # ---- Figure -------------------------------------------------------------------
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
    ax[1].axhline(y.mean(), ls="--", c="k", lw=1, label=f"baseline={y.mean():.3f}")
    ax[1].set(xlabel="Recall", ylabel="Precision", title="OOF Precision-Recall")
    ax[1].legend(loc="upper right")
    top = importance.head(12)[::-1]
    colors = ["crimson" if f in TISSUE else "steelblue" for f in top.index]
    ax[2].barh(top.index, top.values, color=colors)
    ax[2].set(title="Permutation importance (red=tissue feature)", xlabel="mean ROC-AUC drop")
    fig.suptitle(f"Tissue-aware dosage classifier ({best}, {best_combo_name}) "
                f"POS={n_pos} vs TOL={n_neg}", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"\nSaved figure -> {FIG_OUT}")

    # ---- Persist ---------------------------------------------------------------
    import json
    import joblib
    joblib.dump({"pipeline": tuned, "features": FEATURES, "best_family": best,
                "best_combo": best_combo_name}, MODEL_OUT)
    metrics = {
        "n_pos": n_pos, "n_neg": n_neg, "ablation_roc_auc": ablation,
        "best_combo": best_combo_name, "model": best, "best_params": gs.best_params_,
        "cv_roc_auc_tuned": float(gs.best_score_), "oof_roc_auc": float(roc_auc),
        "oof_pr_auc": float(pr_auc), "confusion_matrix": cm.tolist(),
        "permutation_importance": importance.round(5).to_dict(),
    }
    json.dump(metrics, open(METRICS_OUT, "w"), indent=2)
    print(f"Saved model -> {MODEL_OUT}\nSaved metrics -> {METRICS_OUT}")

    # ---- Genome-wide scoring -----------------------------------------------------
    cereb = pd.read_csv(TABLE)
    cereb = cereb[cereb["fit_success"] == True].copy()  # noqa: E712
    sym = symbol_map()
    cereb["symbol"] = cereb["gene"].str.split(".").str[0].map(sym)
    cereb = add_derived(cereb)
    tissue = load_tissue_features()
    genome = cereb.join(tissue, on="symbol", how="inner")
    genome = genome.dropna(subset=FEATURES)
    genome["dosage_score"] = tuned.predict_proba(genome[FEATURES])[:, 1]
    pos = {l.strip() for l in open(POS_FILE) if l.strip()}
    tol_syms = set(pd.read_csv("data/positiveANDnegativeControlGenes.csv", header=None)[0]
                  .dropna().astype(str).str.strip()) - {"geneDUPtol"}
    genome["label"] = np.where(genome["symbol"].isin(pos), "POS",
                              np.where(genome["symbol"].isin(tol_syms), "TOL", ""))
    out = (genome[["gene", "symbol", "dosage_score", "label", "n_obs"]]
          .sort_values("dosage_score", ascending=False).reset_index(drop=True))
    out.to_csv(PRED_OUT, index=False)
    print(f"\nGenome-wide scores -> {PRED_OUT} ({len(out)} genes)")
    print(out.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
