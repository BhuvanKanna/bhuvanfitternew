# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A small toolkit that analyzes per-gene expression distributions from a *C. elegans* natural-isolate dataset (`Supplementary Data 1_csv.csv`). It does two independent things, both per gene:
- fits a **4-parameter Gaussian** + computes **truncation-index** metrics → `fourparam_table.csv`;
- detects **KDE density peaks** (number / location / height / prominence) → `peaks.json`.

## Keeping the GitHub repo updated (important)

This repo (`BhuvanKanna/bhuvanfitternew`, private) is the centralized home for the code **and** the data. **Commit and push every file you change to GitHub as soon as the work in a turn is finished — do this proactively, without being asked.** Do not batch changes for "later." Do not auto-push pre-existing working-tree changes you did not make yourself (e.g. a deletion of a file the user is actively editing) — surface those and confirm first.

Git LFS note: `Supplementary Data 1_csv.csv` (~64 MB) is tracked via **Git LFS** (see `.gitattributes`). Anyone cloning must run `git lfs install` once or the CSV will be a pointer stub. LFS lock-verify can intermittently time out on push; just retry the `git push`.

## Keep this file current

**As you make changes to this project, update `CLAUDE.md` in the same turn** so it always reflects the current architecture, file roles, commands, and the `BhuvanFitter` contract. Treat it as living documentation — when behavior, file layout, or the fit dict changes, adjust the relevant section here (and push it like any other change).

## Single source of truth

`bhuvanfitter.py` holds the analysis library: `_fourparam_gaussian`, the `BhuvanFitter` class, and `gene_peaks` (KDE peak detection). **The notebook and both generator scripts import from it** — never redefine these elsewhere. Make all analysis-logic changes in `bhuvanfitter.py`.

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

A near-identical variant, `generate_fourparam_stats_excluded.py`, fits the same 4-parameter Gaussian but **drops every expression value `<= -0.75`** from each gene's array before fitting (so `n_obs` shrinks and a gene can fall below `MIN_OBS`). It imports `COLUMNS`, `MIN_OBS`, `load_expression`, `git_push`, and `_failed_row` from `generate_fourparam_stats.py` (single source of truth — only `build_table`, the threshold constant `EXCLUDE_AT_OR_BELOW = -0.75`, and output naming differ) and writes/pushes `fourparam_table_excluded.csv`. Same `--limit` / `--no-push` flags.

A parallel pipeline produces the peak dictionary:

```
Supplementary Data 1_csv.csv ──► generate_peaks.py ──► peaks.json
   (load_expression, same loader)   (gene_peaks per col)   (pushed to origin)
```

- `peaks.json` is a nested dict `{gene: {peak_expression_value: {"height":..., "prominence":...}}}`. Peak count per gene = `len()` of the inner dict; a gene with no detectable mode maps to `{}`. **JSON keys are strings**, so the peak expression values are stored as strings — parse to float on load.
- `generate_peaks.py` reuses `load_expression` and `git_push` imported from `generate_fourparam_stats.py` (single definition of each — `git_push` is file-agnostic).

## `BhuvanFitter` contract

- Two fit models are registered in `_FIT_REGISTRY` and dispatched through
  `fit(model, **kwargs)`: **`"fourparam"`** and **`"kde"`**. `active_fits` tracks
  which have run; each fit sets its flag, and `hist(lines=[...])` only draws fits
  that have been run. `**kwargs` are forwarded to the chosen fit method
  (`fourparam` takes none; `kde` takes the `gene_peaks` knobs).
- Histogram is **always 40 bins** (`BhuvanFitter.BINS`).

**`fit("fourparam")`**
- The fit uses `curve_fit(..., method="trf")` with **default linear loss = ordinary least squares**. This is deliberate: it genuinely minimizes the residual sum of squares. (An older Colab version used `loss="soft_l1"`, which is robust but does *not* minimize SSE — it yields different `x0/w` and therefore different metrics. Do not switch back unless intentionally matching that legacy output.)
- `x_max` (constructor arg, defaults to the data max) is the truncation ceiling used by the metrics.
- Returns a 15-key dict:
  `gene, y0, A, x0, w, sumsquarevalue, ti_fourparam_sigma_dist, truncationindex, min, max, right, maxheight, rightheight, n_obs, fit_success`
- `truncationindex` (the renamed height-ratio metric) `== rightheight / maxheight`, where `maxheight = f(peak) − f(min)` and `rightheight = f(x_max) − f(min)`. It returns **NaN** when `maxheight == 0` (degenerate fit whose peak sits at/left of the data minimum — common across the full gene set, so never assume it is finite).
- Metric properties (`truncationindex`, `ti_fourparam_sigma_dist`, `maxheight`, `rightheight`) raise `RuntimeError` until `fit("fourparam")` has been called.

**`fit("kde")`**
- Runs a Gaussian KDE and **reuses the module-level `gene_peaks`** for mode
  detection, so the peaks it reports are identical to `generate_peaks.py` /
  `peaks.json`. Accepts `bw_method`, `min_prominence_frac`, `grid_size`,
  `pad_frac` (same meanings as `gene_peaks`).
- Returns a 6-key dict:
  `gene, n_peaks, peaks, bw_method, n_obs, fit_success`, where `peaks` is the
  `gene_peaks` dict `{value: {"height", "prominence"}}` and `n_peaks == len(peaks)`.
- `fit_success` is `False` (and the cached density is `None`) only when the KDE is
  singular; peak detection still returns `{}` rather than raising. `kde_function(x)`
  evaluates the cached KDE (raises `RuntimeError` if not run / didn't converge).

**Plotting** — `hist(lines=["fourparam", "kde"])` overlays either/both fitted
curves on the 40-bin histogram. The KDE density (which integrates to 1) is scaled
onto the bin-count axis by `n_obs * bin_width`, with `▼` markers at each detected
peak. Unrecognised / not-yet-run / non-converged fits are skipped with a warning.

## `gene_peaks` contract

`gene_peaks(values, min_prominence_frac=0.05, bw_method="silverman", grid_size=1000, pad_frac=0.05, round_to=6)` returns `{peak_expression_value: {"height": <kde density>, "prominence": <prominence>}}` for one gene.

- Detection is on a **Gaussian KDE** (bin-independent), keeping `find_peaks` modes with prominence ≥ `min_prominence_frac` of the max density. Defaults reuse the former `find_density_peaks` settings; tune via the args.
- Returns `{}` for degenerate input (`< 5` finite points, no spread, singular KDE, or no interior mode) — never assume a gene has peaks.
- Peak-value keys are rounded to `round_to` decimals.

## Commands

There is no build system, linter, or test suite. Development is `python` scripts + the notebook. Requires `numpy`, `pandas`, `scipy`, `matplotlib`.

```bash
# Regenerate the full table from the CSV and push it to the repo (~25.8k genes)
python generate_fourparam_stats.py

# Quick sanity check: first 50 genes, write the CSV but DON'T push
python generate_fourparam_stats.py --limit 50 --no-push

# Generate the full table locally without pushing
python generate_fourparam_stats.py --no-push

# Peak dictionary — same flags (all genes + push / sanity subset / local only)
python generate_peaks.py
python generate_peaks.py --limit 50 --no-push
python generate_peaks.py --no-push
```

Both generators take `--limit N` and `--no-push`. After a `--limit` run the output file (`fourparam_table.csv` or `peaks.json`) holds only those genes; restore the committed full version with `git checkout -- <file>` before pushing anything.

The notebook (`newbhuvanfitter.ipynb`) is just `from bhuvanfitter import BhuvanFitter, gene_peaks` plus a synthetic single-gene example — use it for interactive inspection of one gene's `fit("fourparam")` / `fit("kde")` and the `hist(lines=["fourparam", "kde"])` overlay.

## Other data files

`Supplementary Data 1 trunc 20250702.xlsx` (gene-name ↔ identifier mapping) and `genes_of_interest.json` (curated gene sets) are inputs for downstream gene-of-interest analysis; they are not consumed by `generate_fourparam_stats.py`.
