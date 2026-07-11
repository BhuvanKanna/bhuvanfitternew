#!/usr/bin/env python
"""Combine level + pure-shape features with a separate regularization strength
for each group, and compare against the uniform-penalty baselines.

Plain sklearn LogisticRegression applies one C to every (standardized) feature
equally. To let LEVEL and PURE_SHAPE get different effective penalties without
a custom loss, this scales the shape columns by a tunable factor `scale` AFTER
standardization, then fits a single global C. Scaling a feature down forces the
optimizer to use a larger raw coefficient for the same effect, and that raw
coefficient is what the L2 penalty sees -- so scale < 1 tightens the shape
group's effective penalty relative to level, scale > 1 relaxes it, and scale=1
recovers the ordinary uniform-penalty combined model. Grid-searching (C, scale)
jointly explores the full space of relative level/shape regularization.

Run: python train_grouped_regularization.py
"""
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from ablation_dosage_features import LEVEL, PURE_SHAPE
from train_dosage_classifier import load_labelled

FEATURES = LEVEL + PURE_SHAPE
GROUP_MASK = np.array([False] * len(LEVEL) + [True] * len(PURE_SHAPE))


class GroupScaler(BaseEstimator, TransformerMixin):
    """Multiply the columns flagged by `mask` by `scale`; leave the rest as-is."""

    def __init__(self, mask=None, scale=1.0):
        self.mask = mask
        self.scale = scale

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X = np.asarray(X, dtype=float).copy()
        if self.mask is not None:
            X[:, self.mask] *= self.scale
        return X


def make_grouped_pipeline():
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("groupscale", GroupScaler(mask=GROUP_MASK)),
        ("clf", LogisticRegression(max_iter=5000, class_weight="balanced")),
    ])


def main():
    lab = load_labelled()
    X, y = lab[FEATURES], lab["label"].values
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=0)
    print(f"POS={int(y.sum())}  TOL={int((y==0).sum())}  prevalence={y.mean():.3f}")
    print(f"Features: {len(LEVEL)} level + {len(PURE_SHAPE)} pure-shape = {len(FEATURES)}\n")

    # ---- Uniform-penalty baseline (scale=1, i.e. today's ordinary combined model) ----
    uni_grid = {"clf__C": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]}
    uni = GridSearchCV(make_grouped_pipeline(), uni_grid, scoring="roc_auc", cv=cv, n_jobs=-1)
    uni.fit(X, y)
    uni_pa = cross_validate(uni.best_estimator_, X, y, cv=cv, scoring="average_precision",
                            n_jobs=-1)["test_score"].mean()
    print(f"Uniform-penalty combined (scale=1): best C={uni.best_params_['clf__C']}  "
          f"ROC-AUC={uni.best_score_:.3f}  PR-AUC={uni_pa:.3f}")

    # ---- Joint (C, group scale) search: separate effective penalty per group ----
    grid = {
        "clf__C": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
        "groupscale__scale": [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0],
    }
    gs = GridSearchCV(make_grouped_pipeline(), grid, scoring="roc_auc", cv=cv, n_jobs=-1)
    gs.fit(X, y)
    best_pa = cross_validate(gs.best_estimator_, X, y, cv=cv, scoring="average_precision",
                             n_jobs=-1)["test_score"].mean()
    C, scale = gs.best_params_["clf__C"], gs.best_params_["groupscale__scale"]
    rel_penalty = 1.0 / (scale ** 2)  # shape penalty relative to level (1.0 = equal)
    print(f"\nGrouped-penalty combined: best C={C}  shape_scale={scale}  "
          f"(shape penalized {rel_penalty:.3g}x relative to level)")
    print(f"ROC-AUC={gs.best_score_:.3f}  PR-AUC={best_pa:.3f}")

    print(f"\nFor reference: level_only=0.818, mean_only=0.824, "
          f"pure_shape(own best C)=0.775 (see ablation_dosage_features.py)")
    delta = gs.best_score_ - uni.best_score_
    print(f"\nGrouped vs uniform-penalty combined: {'+' if delta >= 0 else ''}{delta:.3f} ROC-AUC")

    # ---- Full grid, sorted, for transparency ----
    import pandas as pd
    res = pd.DataFrame(gs.cv_results_)[["param_clf__C", "param_groupscale__scale",
                                        "mean_test_score", "std_test_score"]]
    res = res.sort_values("mean_test_score", ascending=False).head(10)
    print("\nTop 10 (C, shape_scale) combinations by CV ROC-AUC:")
    print(res.rename(columns={"param_clf__C": "C", "param_groupscale__scale": "shape_scale",
                              "mean_test_score": "roc_auc", "std_test_score": "std"})
          .to_string(index=False))


if __name__ == "__main__":
    main()
