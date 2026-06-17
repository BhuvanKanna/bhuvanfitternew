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

**Do this automatically and proactively for every significant change — the user should not have to ask.** A significant change is anything that makes a section here stale: behavior/logic changes, new or removed files, fit-dict keys, notebook cells, new dependencies, or command/flag changes. After any non-trivial edit, check whether the relevant section is now out of date and, if so, fix it in the same commit/push.

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

A near-identical variant, `generate_fourparam_stats_excluded.py`, fits the same 4-parameter Gaussian but **drops every expression value `<= EXCLUDE_AT_OR_BELOW`** (currently `-1`) from each gene's array before fitting (so `n_obs` shrinks and a gene can fall below `MIN_OBS`). It imports `COLUMNS`, `MIN_OBS`, `load_expression`, `git_push`, and `_failed_row` from `generate_fourparam_stats.py` (single source of truth — only `build_table`, the threshold constant `EXCLUDE_AT_OR_BELOW`, and output naming differ). **The output filename encodes the threshold** — `OUTPUT_CSV = HERE / f"fourparam_table_excluded_at_or_below_{EXCLUDE_AT_OR_BELOW}.csv"` (and the commit message is likewise threshold-derived) — so changing `EXCLUDE_AT_OR_BELOW` writes/pushes a *separate*, self-labeled spreadsheet (e.g. `fourparam_table_excluded_at_or_below_-1.csv`) instead of overwriting a differently-thresholded one. Same `--limit` / `--no-push` flags. The earlier `-0.75` run is preserved as `fourparam_table_excluded_at_or_below_-0.75.csv` (regenerate by setting `EXCLUDE_AT_OR_BELOW = -0.75`).

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
- `truncationindex` (the renamed height-ratio metric) `== rightheight / maxheight`, where the baseline subtracted from both is the **fitted curve's minimum over the histogram interval** (`_curve_baseline()`, i.e. `min(f)` on a 600-pt grid over `[hist_edges[0], hist_edges[-1]]`) — **not** `f(data min)`. So `maxheight = max(f) − min(f)` and `rightheight = f(x_max) − min(f)`. Because the baseline is the curve's true interval-minimum, the ratio is **bounded to [0, 1]** (0 = ceiling at the curve minimum, 1 = ceiling at the peak). It returns **NaN** only when `maxheight == 0` (the curve is flat over the interval). (An earlier version used `f(data min)` as the baseline, which left the ratio unbounded — values ranged to ±1e14 — so don't reintroduce that.)
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

There is no build system, linter, or test suite. Development is `python` scripts + the notebook. Requires `numpy`, `pandas`, `scipy`, `matplotlib` (plus `openpyxl` if you read `Supplementary Data 1 trunc 20250702.xlsx` via `pd.read_excel`).

```bash
# Regenerate the full table from the CSV and push it to the repo (~25.8k genes)
python generate_fourparam_stats.py

# Quick sanity check: first 50 genes, write the CSV but DON'T push
python generate_fourparam_stats.py --limit 50 --no-push

# Generate the full table locally without pushing
python generate_fourparam_stats.py --no-push

# Excluded variant: default threshold -1, or override per run (filename encodes it)
python generate_fourparam_stats_excluded.py                    # -> ..._at_or_below_-1.csv (push)
python generate_fourparam_stats_excluded.py --threshold -0.75  # -> ..._at_or_below_-0.75.csv
python generate_fourparam_stats_excluded.py --limit 50 --no-push

# Peak dictionary — same flags (all genes + push / sanity subset / local only)
python generate_peaks.py
python generate_peaks.py --limit 50 --no-push
python generate_peaks.py --no-push
```

All generators take `--limit N` and `--no-push`; `generate_fourparam_stats_excluded.py` additionally takes `--threshold T` (default `-1`) which sets the exclusion cutoff **and** the output filename (`fourparam_table_excluded_at_or_below_<T>.csv`). After a `--limit` run the output file holds only those genes; restore the committed full version with `git checkout -- <file>` before pushing anything. To regenerate several thresholds at once, run each with `--no-push` (distinct output files, no git race) and make one commit.

The notebook (`newbhuvanfitter.ipynb`) is the interactive scratch space (`from bhuvanfitter import BhuvanFitter, gene_peaks`). It does three things:
- **Per-gene inspection** — build a `BhuvanFitter` on one column of `master` (the transposed `Supplementary Data 1_csv.csv`) and view its `fit("fourparam")` / `fit("kde")` and the `hist(lines=["fourparam", "kde"])` overlay.
- **Cross-gene distribution plots** — it loads `fourparam_table_excluded_at_or_below_-1.csv` into `fourparam_df` and defines two helpers for histogramming any column across genes: `select(param)` applies the **single shared filter** (`fit_success == True` and `0 < truncationindex < 1`, NaNs dropped) and `plot_param_hist(param, *, color, bins, log)` draws the histogram (returning the data) and prints a filtering funnel — gene counts at each stage (master total → fourparam_df total → fit_success → NaN drop → truncationindex > 0 → < 1). Change the filter in one place (`select`) and both the sumsquarevalue and truncationindex plots follow.
- **Per-strain expression summaries (filtered)** — the last three cells restrict `master` to the genes passing the **same filter** (`fit_success == True` and `0 < truncationindex < 1`, ~11k genes) via `filtered_master = master[fourparam_df.loc[filter_mask, 'gene']]`, then histogram the per-column **mean** (`strain_avgs = filtered_master.mean(axis=0)`) and **standard deviation** (`strain_sds = filtered_master.std(axis=0)`) across the 207 isolates, plus a **mean-vs-SD scatter** (`strain_avgs` x, `strain_sds` y; aligned by gene index). The SD and scatter cells reuse `filtered_master` / `strain_avgs` / `strain_sds` defined in the mean cell, so run them in order.

## Other data files

`Supplementary Data 1 trunc 20250702.xlsx` (gene-name ↔ identifier mapping) and `genes_of_interest.json` (curated gene sets) are inputs for downstream gene-of-interest analysis; they are not consumed by `generate_fourparam_stats.py`.
