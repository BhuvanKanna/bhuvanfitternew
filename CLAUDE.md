# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A small toolkit that fits a **4-parameter Gaussian** to per-gene expression histograms from a *C. elegans* natural-isolate dataset and computes **truncation-index** metrics. It produces `fourparam_table.csv` (one row per gene) from `Supplementary Data 1_csv.csv`.

## Keeping the GitHub repo updated (important)

This repo (`BhuvanKanna/bhuvanfitternew`, private) is the centralized home for the code **and** the data. **Commit and push every file you change to GitHub as soon as the work in a turn is finished — do this proactively, without being asked.** Do not batch changes for "later." Do not auto-push pre-existing working-tree changes you did not make yourself (e.g. a deletion of a file the user is actively editing) — surface those and confirm first.

Git LFS note: `Supplementary Data 1_csv.csv` (~64 MB) is tracked via **Git LFS** (see `.gitattributes`). Anyone cloning must run `git lfs install` once or the CSV will be a pointer stub. LFS lock-verify can intermittently time out on push; just retry the `git push`.

## Single source of truth

`bhuvanfitter.py` holds the `_fourparam_gaussian` function and the `BhuvanFitter` class. **Both `newbhuvanfitter.ipynb` and `generate_fourparam_stats.py` import from it** — never redefine the class in either of those. Make all fitting-logic changes in `bhuvanfitter.py`.

## The pipeline (requires reading 3 files to see end-to-end)

```
Supplementary Data 1_csv.csv          bhuvanfitter.py
  (strain col + 207 isolate cols,       (BhuvanFitter.fit)
   25,849 gene rows)                           │
        │  set_index('strain').T               │
        ▼  → 207 strains × 25,849 genes         ▼
  generate_fourparam_stats.py ── per gene ──► fourparam_table.csv
        │   (build_table loops df.columns)      (25,849 rows, 15 cols)
        └── commits + pushes the CSV to origin
```

- Genes with `< 10` finite observations, or whose `curve_fit` fails to converge, are written as a row with `fit_success=False` and NaN metrics rather than skipped or crashing.
- `generate_fourparam_stats.py` has a module-level `COLUMNS` list that **must stay identical to the keys `BhuvanFitter.fit("fourparam")` returns** (same names, same order). If you add/rename a key in the fit dict, update `COLUMNS`.

## `BhuvanFitter` contract

- Histogram is **always 40 bins** (`BhuvanFitter.BINS`).
- The fit uses `curve_fit(..., method="trf")` with **default linear loss = ordinary least squares**. This is deliberate: it genuinely minimizes the residual sum of squares. (An older Colab version used `loss="soft_l1"`, which is robust but does *not* minimize SSE — it yields different `x0/w` and therefore different metrics. Do not switch back unless intentionally matching that legacy output.)
- `x_max` (constructor arg, defaults to the data max) is the truncation ceiling used by the metrics.
- `fit("fourparam")` returns a 15-key dict:
  `gene, y0, A, x0, w, sumsquarevalue, ti_fourparam_sigma_dist, truncationindex, min, max, right, maxheight, rightheight, n_obs, fit_success`
- `truncationindex` (the renamed height-ratio metric) `== rightheight / maxheight`, where `maxheight = f(peak) − f(min)` and `rightheight = f(x_max) − f(min)`. It returns **NaN** when `maxheight == 0` (degenerate fit whose peak sits at/left of the data minimum — common across the full gene set, so never assume it is finite).
- Metric properties (`truncationindex`, `ti_fourparam_sigma_dist`, `maxheight`, `rightheight`) raise `RuntimeError` until `fit` has been called.

## Commands

There is no build system, linter, or test suite. Development is `python` scripts + the notebook. Requires `numpy`, `pandas`, `scipy`, `matplotlib`.

```bash
# Regenerate the full table from the CSV and push it to the repo (~25.8k genes)
python generate_fourparam_stats.py

# Quick sanity check: first 50 genes, write the CSV but DON'T push
python generate_fourparam_stats.py --limit 50 --no-push

# Generate the full table locally without pushing
python generate_fourparam_stats.py --no-push
```

After a `--limit` run, `fourparam_table.csv` is truncated to that many rows; restore the committed full table with `git checkout -- fourparam_table.csv` before pushing anything.

The notebook (`newbhuvanfitter.ipynb`) is just `from bhuvanfitter import BhuvanFitter` plus a synthetic single-gene example — use it for interactive inspection of one gene's fit and `hist()` plot.

## Other data files

`Supplementary Data 1 trunc 20250702.xlsx` (gene-name ↔ identifier mapping) and `genes_of_interest.json` (curated gene sets) are inputs for downstream gene-of-interest analysis; they are not consumed by `generate_fourparam_stats.py`.
