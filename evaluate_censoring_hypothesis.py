#!/usr/bin/env python
"""Directly test the grant's "Censoring Hypothesis" with a level-matched evaluation.

Claim (grant.pdf / acfrog.md): a gene's natural expression-level distribution
across a healthy population should be Gaussian. If a gene has a deleterious
overexpression (OE) threshold, individuals who'd express past it are
effectively censored from the "healthy" sample, truncating the right tail
(an expression "ceiling"). So OE-sensitive genes (POS) should show a
truncated/capped distribution shape; OE-tolerant genes (TOL) an uncapped one.
Tested in the proposal in two species -- wild C. elegans and human GTEx tissue.

The repo's existing dosage classifier already shows shape features carry some
signal in isolation, but nothing on top of expression *level* -- a linear
(regularization-based) test that can still let level leak through nonlinearly.
This script tests the claim more directly: match each POS gene to its nearest
TOL genes BY EXPRESSION LEVEL, then ask whether shape/truncation still
separates POS from TOL within that level-matched set, in both species.

Run: python evaluate_censoring_hypothesis.py
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.stats import kurtosis, mannwhitneyu, skew
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict, cross_validate
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore", category=UserWarning)

from ablation_dosage_features import PURE_SHAPE
from generate_fourparam_stats import load_expression
from regenerate_acfrog_figures import EXCLUDE_AT_OR_BELOW
from regenerate_grant_figures import (
    build_worm_groups_names,
    human_to_worm,
    load_pos,
    load_tol,
    _worm_table,
)
from train_dosage_classifier import load_labelled as load_labelled_human
from train_dosage_classifier import make_pipeline

K_NEIGHBORS = 5
FIG_OUT = "censoring_hypothesis_eval.png"
SUMMARY_OUT = "censoring_hypothesis_summary.csv"


# ---------------------------------------------------------------------------
# Worm labelled table (mirrors train_dosage_classifier.load_labelled, but for
# worm: ortholog-bridged POS/TOL, canonical transcript per gene, and a
# no-refit skew/kurt/mean/std pass since the worm table lacks those columns).
# ---------------------------------------------------------------------------
def load_labelled_worm():
    pos_h, tol_h = load_pos(), load_tol()
    orth = human_to_worm(pos_h | tol_h)
    worm_pos_names = {w for h in pos_h for w in orth.get(h, [])}
    worm_tol_names = {w for h in tol_h for w in orth.get(h, [])} - worm_pos_names

    ids_pos, _ = build_worm_groups_names(worm_pos_names)
    ids_tol, _ = build_worm_groups_names(worm_tol_names)

    tab = _worm_table()
    tab = tab[tab["fit_success"] == True]  # noqa: E712

    def canonical(ids, label):
        sub = tab.loc[[i for i in ids if i in tab.index]].copy()
        sub = sub.sort_values("n_obs", ascending=False)
        sub = sub[~sub.index.duplicated(keep="first")]  # dedup any repeated ids
        sub["label"] = label
        return sub

    lab = pd.concat([canonical(ids_pos, 1), canonical(ids_tol, 0)])
    lab = lab.reset_index().rename(columns={"index": "gene"})

    # No-refit skew/kurt/mean/std: pull raw values from worm.csv for just
    # these genes (mirrors the cerebellum notebook's build_shape_features).
    master = load_expression("worm.csv")
    means, stds, sks, kus = [], [], [], []
    for g in lab["gene"]:
        v = master[g].astype(float).values
        v = v[np.isfinite(v)]
        v = v[v > EXCLUDE_AT_OR_BELOW]
        if len(v) < 3:
            means.append(np.nan); stds.append(np.nan)
            sks.append(np.nan); kus.append(np.nan)
            continue
        means.append(v.mean()); stds.append(v.std(ddof=1))
        sks.append(skew(v)); kus.append(kurtosis(v))
    lab["mean"], lab["std"], lab["skew"], lab["kurt"] = means, stds, sks, kus
    lab["cv"] = lab["std"] / lab["mean"]
    lab["w_x0"] = lab["w"] / lab["x0"]
    return lab.replace([np.inf, -np.inf], np.nan)


# ---------------------------------------------------------------------------
# Level-matched evaluation
# ---------------------------------------------------------------------------
def match_by_level(lab, k=K_NEIGHBORS):
    """Match each POS gene to its k nearest TOL genes by `mean` (with
    replacement). Returns the matched sub-dataframe (all POS + unique matched TOL)."""
    lab = lab.dropna(subset=["mean"] + PURE_SHAPE).reset_index(drop=True)
    pos = lab[lab["label"] == 1]
    tol = lab[lab["label"] == 0]
    nn = NearestNeighbors(n_neighbors=min(k, len(tol))).fit(tol[["mean"]])
    _, idx = nn.kneighbors(pos[["mean"]])
    matched_tol_idx = sorted(set(idx.ravel()))
    matched_tol = tol.iloc[matched_tol_idx]
    return pd.concat([pos, matched_tol]).reset_index(drop=True), lab


def shape_auc(df):
    """CV ROC-AUC/PR-AUC of a PURE_SHAPE-only logistic regression on df."""
    X, y = df[PURE_SHAPE], df["label"].values
    if y.sum() < 2 or (y == 0).sum() < 2:
        return np.nan, np.nan, None
    n_splits = min(5, int(y.sum()), int((y == 0).sum()))
    if n_splits < 2:
        return np.nan, np.nan, None
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=10, random_state=0)
    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0),
                         scale=True)
    res = cross_validate(pipe, X, y, cv=cv, scoring=["roc_auc", "average_precision"], n_jobs=-1)
    oof_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=1, random_state=0)
    oof = cross_val_predict(pipe, X, y, cv=oof_cv, method="predict_proba", n_jobs=-1)[:, 1]
    return res["test_roc_auc"].mean(), res["test_average_precision"].mean(), (y, oof)


def evaluate(species, lab):
    print(f"\n{'='*70}\n{species.upper()}\n{'='*70}")
    print(f"Labelled genes: POS={int((lab['label']==1).sum())}  "
          f"TOL={int((lab['label']==0).sum())}")

    matched, lab_clean = match_by_level(lab)
    pos = matched[matched["label"] == 1]
    tol_matched = matched[matched["label"] == 0]
    tol_all = lab_clean[lab_clean["label"] == 0]

    print(f"Matched set: POS={len(pos)}  matched-TOL={len(tol_matched)} "
          f"(from {len(tol_all)} total TOL, k={K_NEIGHBORS})")
    print(f"Median mean-expression -- POS: {pos['mean'].median():.3f}  "
          f"all-TOL: {tol_all['mean'].median():.3f}  "
          f"matched-TOL: {tol_matched['mean'].median():.3f}  "
          f"(matched should sit close to POS)")

    # Directional prediction from the Censoring Hypothesis applies cleanly only
    # to the two features this repo's BhuvanFitter contract defines as
    # truncation-strength metrics: higher truncationindex = more capped (POS
    # should be greater); lower ti_fourparam_sigma_dist = ceiling closer to the
    # peak, i.e. more capped (POS should be LESS). skew/kurt/cv/w_x0 have no
    # unambiguous a-priori direction from the hypothesis, so those get a
    # two-sided test instead of an arbitrary "greater".
    DIRECTION = {"truncationindex": "greater", "ti_fourparam_sigma_dist": "less"}
    print(f"\n{'feature':24s} {'alt.':>10s} {'MWU p (matched)':>18s}")
    for feat in PURE_SHAPE:
        a, b = pos[feat].dropna(), tol_matched[feat].dropna()
        alt = DIRECTION.get(feat, "two-sided")
        if len(a) >= 3 and len(b) >= 3:
            p = mannwhitneyu(a, b, alternative=alt).pvalue
        else:
            p = np.nan
        print(f"{feat:24s} {alt:>10s} {p:18.4f}")

    unmatched_auc, unmatched_pr, _ = shape_auc(lab_clean)
    matched_auc, matched_pr, matched_oof = shape_auc(matched)
    print(f"\nShape-only (PURE_SHAPE) ROC-AUC/PR-AUC:")
    print(f"  Unmatched (full population): ROC-AUC={unmatched_auc:.3f}  PR-AUC={unmatched_pr:.3f}")
    print(f"  Level-matched              : ROC-AUC={matched_auc:.3f}  PR-AUC={matched_pr:.3f}")

    ti_p = mannwhitneyu(pos["truncationindex"].dropna(), tol_matched["truncationindex"].dropna(),
                        alternative="greater").pvalue
    sigma_p = mannwhitneyu(pos["ti_fourparam_sigma_dist"].dropna(),
                           tol_matched["ti_fourparam_sigma_dist"].dropna(),
                           alternative="less").pvalue
    verdict = ("SURVIVES -- shape/truncation still separates POS from level-matched TOL"
              if matched_auc >= 0.65 and (ti_p < 0.05 or sigma_p < 0.05) else
              "DOES NOT SURVIVE -- shape signal collapses once expression level is controlled for")
    print(f"truncationindex p={ti_p:.4f} (greater)   ti_fourparam_sigma_dist p={sigma_p:.4f} (less)")
    print(f"\nVerdict ({species}): {verdict}")

    return {
        "species": species, "n_pos": len(pos), "n_tol_all": len(tol_all),
        "n_tol_matched": len(tol_matched),
        "median_mean_pos": pos["mean"].median(),
        "median_mean_tol_all": tol_all["mean"].median(),
        "median_mean_tol_matched": tol_matched["mean"].median(),
        "median_ti_pos": pos["truncationindex"].median(),
        "median_ti_tol_matched": tol_matched["truncationindex"].median(),
        "mwu_p_ti_matched": ti_p,
        "unmatched_shape_roc_auc": unmatched_auc, "unmatched_shape_pr_auc": unmatched_pr,
        "matched_shape_roc_auc": matched_auc, "matched_shape_pr_auc": matched_pr,
        "verdict": verdict,
    }, matched_oof


def main():
    lab_h = load_labelled_human()
    lab_w = load_labelled_worm()

    rows = []
    oofs = {}
    for species, lab in [("human_cerebellum", lab_h), ("worm", lab_w)]:
        row, oof = evaluate(species, lab)
        rows.append(row)
        oofs[species] = oof

    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_OUT, index=False)
    print(f"\nSummary -> {SUMMARY_OUT}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for row_i, (species, lab) in enumerate([("human_cerebellum", lab_h), ("worm", lab_w)]):
        matched, lab_clean = match_by_level(lab)
        pos, tol_all = lab_clean[lab_clean["label"] == 1], lab_clean[lab_clean["label"] == 0]
        tol_matched = matched[matched["label"] == 0]

        ax = axes[row_i, 0]
        bins = np.linspace(lab_clean["mean"].min(), lab_clean["mean"].max(), 30)
        ax.hist(tol_all["mean"], bins=bins, alpha=0.35, color="gray", label="TOL (all)", density=True)
        ax.hist(tol_matched["mean"], bins=bins, alpha=0.55, color="steelblue",
                label="TOL (matched)", density=True)
        ax.hist(pos["mean"], bins=bins, alpha=0.55, color="crimson", label="POS", density=True)
        ax.set(title=f"{species}: mean expression, POS vs TOL", xlabel="mean expression")
        ax.legend(fontsize=8)

        ax2 = axes[row_i, 1]
        y, oof = oofs[species]
        if oof is not None:
            fpr, tpr, _ = roc_curve(y, oof)
            auc = roc_auc_score(y, oof)
            ax2.plot(fpr, tpr, lw=2, color="purple", label=f"matched shape AUC={auc:.3f}")
        ax2.plot([0, 1], [0, 1], "k--", lw=1)
        ax2.set(title=f"{species}: level-matched shape-only ROC", xlabel="FPR", ylabel="TPR")
        ax2.legend(fontsize=8, loc="lower right")

    fig.suptitle("Censoring Hypothesis: level-matched shape/truncation evaluation", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"Saved figure -> {FIG_OUT}")


if __name__ == "__main__":
    main()
