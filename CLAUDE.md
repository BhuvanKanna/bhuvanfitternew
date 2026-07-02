# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A small toolkit that analyzes per-gene expression distributions from a *C. elegans* natural-isolate dataset (`worm.csv`). It does two independent things, both per gene:
- fits a **4-parameter Gaussian** + computes **truncation-index** metrics → `worm_fourparam_table.csv`;
- detects **KDE density peaks** (number / location / height / prominence) → `peaks.json`.

## Keeping the GitHub repo updated (important)

This repo (`BhuvanKanna/bhuvanfitternew`, private) is the centralized home for the code **and** the data. **Commit and push every file you change to GitHub as soon as the work in a turn is finished — do this proactively, without being asked.** Do not batch changes for "later." Do not auto-push pre-existing working-tree changes you did not make yourself (e.g. a deletion of a file the user is actively editing) — surface those and confirm first.

Git LFS note: the large source matrices — `worm.csv` (~64 MB, worm; originally `Supplementary Data 1_csv.csv`) and `cerebellumlog2.csv` (~203 MB, GTEx cerebellum) — are tracked via **Git LFS** (see `.gitattributes`). Anyone cloning must run `git lfs install` once or these CSVs will be pointer stubs. LFS lock-verify can intermittently time out on push; just retry the `git push`. (The generated per-gene tables, e.g. `cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv`, are committed as normal files.)

## Keep this file current

**As you make changes to this project, update `CLAUDE.md` in the same turn** so it always reflects the current architecture, file roles, commands, and the `BhuvanFitter` contract. Treat it as living documentation — when behavior, file layout, or the fit dict changes, adjust the relevant section here (and push it like any other change).

**Do this automatically and proactively for every significant change — the user should not have to ask.** A significant change is anything that makes a section here stale: behavior/logic changes, new or removed files, fit-dict keys, notebook cells, new dependencies, or command/flag changes. After any non-trivial edit, check whether the relevant section is now out of date and, if so, fix it in the same commit/push.

## Single source of truth

`bhuvanfitter.py` holds the analysis library: `_fourparam_gaussian`, the `BhuvanFitter` class, and `gene_peaks` (KDE peak detection). **The notebook and both generator scripts import from it** — never redefine these elsewhere. Make all analysis-logic changes in `bhuvanfitter.py`.

## The pipeline (requires reading 3 files to see end-to-end)

```
worm.csv                              bhuvanfitter.py
  (strain col + 207 isolate cols,       (BhuvanFitter.fit)
   25,849 gene rows)                           │
        │  set_index('strain').T               │
        ▼  → 207 strains × 25,849 genes         ▼
  generate_fourparam_stats.py ── per gene ──► worm_fourparam_table.csv
        │   (build_table loops df.columns)      (25,849 rows, 17 cols)
        └── commits + pushes the CSV to origin
```

- Genes with `< 10` finite observations, or whose `curve_fit` fails to converge, are written as a row with `fit_success=False` and NaN metrics rather than skipped or crashing.
- `generate_fourparam_stats.py` has a module-level `COLUMNS` list that **must stay identical to the keys `BhuvanFitter.fit("fourparam")` returns** (same names, same order). If you add/rename a key in the fit dict, update `COLUMNS`.
- `load_expression(csv_path, id_col="strain", drop_cols=())` is the shared loader (genes as rows → transposed to samples × genes). It is **dataset-agnostic**: `id_col` names the gene-identifier column (the worm file's first column is mislabeled `strain` but holds gene names; GTEx uses `Name`) and `drop_cols` discards extra non-sample index columns (e.g. GTEx's `Description`). Defaults reproduce the original worm behaviour exactly.

A near-identical variant, `generate_fourparam_stats_excluded.py`, fits the same 4-parameter Gaussian but **drops every expression value `<= EXCLUDE_AT_OR_BELOW`** (default `-1`) from each gene's array before fitting (so `n_obs` shrinks and a gene can fall below `MIN_OBS`). It imports `COLUMNS`, `MIN_OBS`, `load_expression`, `git_push`, and `_failed_row` from `generate_fourparam_stats.py` (single source of truth — only `build_table` / the parallel worker, the threshold constant `EXCLUDE_AT_OR_BELOW`, and output naming differ). **The output filename encodes the input dataset and the threshold** via `output_csv_for(threshold, input_path)` — `<input stem>_fourparam_table_excluded_at_or_below_<threshold>.csv` — so datasets never collide: worm → `worm_fourparam_table_excluded_at_or_below_-1.csv`, GTEx → `cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv`. The earlier worm `-0.75` run is preserved as `worm_fourparam_table_excluded_at_or_below_-0.75.csv`.

  - **Flags:** `--threshold T` (exclusion cutoff, default `-1`), `--limit N`, `--no-push`, plus `--input PATH` / `--id-col COL` / `--drop-col COL` (repeatable) for targeting other datasets through the generalized `load_expression`, and **`--jobs N`** (parallel worker processes; genes are independent so this scales near-linearly) and **`--max-nfev N`** (curve_fit cap, forwarded to `fit("fourparam")`). Parallelism uses a `multiprocessing.Pool` whose workers cache the threshold / `max_nfev` via `_init_worker`; the module-level `_fit_one((gene, values))` worker is picklable and `imap` preserves input order, so a parallel run is **bit-identical** to `--jobs 1`. The GTEx cerebellum table was generated with `--input cerebellumlog2.csv --id-col Name --drop-col Description --threshold -1 --jobs 11 --max-nfev 2000`.

A parallel pipeline produces the peak dictionary:

```
worm.csv ──► generate_peaks.py ──► peaks.json
   (load_expression, same loader)   (gene_peaks per col)   (pushed to origin)
```

- `peaks.json` is a nested dict `{gene: {peak_expression_value: {"height":..., "prominence":...}}}`. Peak count per gene = `len()` of the inner dict; a gene with no detectable mode maps to `{}`. **JSON keys are strings**, so the peak expression values are stored as strings — parse to float on load.
- `generate_peaks.py` reuses `load_expression` and `git_push` imported from `generate_fourparam_stats.py` (single definition of each — `git_push` is file-agnostic).

## `BhuvanFitter` contract

- Two fit models are registered in `_FIT_REGISTRY` and dispatched through
  `fit(model, **kwargs)`: **`"fourparam"`** and **`"kde"`**. `active_fits` tracks
  which have run; each fit sets its flag, and `hist(lines=[...])` only draws fits
  that have been run. `**kwargs` are forwarded to the chosen fit method
  (`fourparam` takes `max_nfev`; `kde` takes the `gene_peaks` knobs).
- Histogram is **always 40 bins** (`BhuvanFitter.BINS`).

**`fit("fourparam")`**
- The fit uses `curve_fit(..., method="trf")` with **default linear loss = ordinary least squares**. This is deliberate: it genuinely minimizes the residual sum of squares. (An older Colab version used `loss="soft_l1"`, which is robust but does *not* minimize SSE — it yields different `x0/w` and therefore different metrics. Do not switch back unless intentionally matching that legacy output.)
- `x_max` (constructor arg, defaults to the data max) is the truncation ceiling used by the metrics.
- `fit("fourparam", max_nfev=...)` caps the curve_fit function evaluations (default `10_000`). Genuine fits converge well under it; lowering it (e.g. `2000`) mainly cuts wasted effort on non-converging genes — which otherwise burn the whole budget before raising — at the small risk of flipping a slow-converging gene to `fit_success=False`. The generators expose it via `--max-nfev`.
- Returns a 17-key dict:
  `gene, y0, A, x0, w, sumsquarevalue, ti_fourparam_sigma_dist, truncationindex, min, max, mean, std, right, maxheight, rightheight, n_obs, fit_success`
- `mean` and `std` are pure per-gene summary statistics of the finite values used for the fit (`mean = data.mean()`, `std` = **sample** standard deviation, `ddof=1`; NaN when `n_obs <= 1`). For the excluded variant they reflect the post-exclusion array. They sit alongside `min`/`max`/`n_obs` and do not depend on whether the fit converged (but a `fit_success=False` row still writes them as NaN, like every other metric).
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

There is no build system, linter, or test suite. Development is `python` scripts + the notebook. Requires `numpy`, `pandas`, `scipy`, `matplotlib` (plus `openpyxl` if you read `Supplementary Data 1 trunc 20250702.xlsx` via `pd.read_excel`, and `scikit-learn` for `cerebellumbhuvanfitter.ipynb` section 7's pTriplo-prediction model).

```bash
# Regenerate the full table from the CSV and push it to the repo (~25.8k genes)
python generate_fourparam_stats.py

# Quick sanity check: first 50 genes, write the CSV but DON'T push
python generate_fourparam_stats.py --limit 50 --no-push

# Generate the full table locally without pushing
python generate_fourparam_stats.py --no-push

# Excluded variant: default threshold -1, or override per run (filename encodes it)
python generate_fourparam_stats_excluded.py                    # -> worm_..._at_or_below_-1.csv (push)
python generate_fourparam_stats_excluded.py --threshold -0.75  # -> worm_..._at_or_below_-0.75.csv
python generate_fourparam_stats_excluded.py --limit 50 --no-push

# Excluded variant on another dataset (GTEx cerebellum), parallel + lower curve_fit cap
python generate_fourparam_stats_excluded.py --input cerebellumlog2.csv \
    --id-col Name --drop-col Description --threshold -1 --jobs 11 --max-nfev 2000 --no-push
    # -> cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv

# Peak dictionary — same flags (all genes + push / sanity subset / local only)
python generate_peaks.py
python generate_peaks.py --limit 50 --no-push
python generate_peaks.py --no-push
```

All generators take `--limit N` and `--no-push`; `generate_fourparam_stats_excluded.py` additionally takes `--threshold T` (default `-1`, sets the exclusion cutoff **and** the output filename), `--input` / `--id-col` / `--drop-col` (target another dataset), and `--jobs N` / `--max-nfev N` (parallelism and the curve_fit cap). After a `--limit` run the output file holds only those genes; restore the committed full version with `git checkout -- <file>` before pushing anything. To regenerate several thresholds at once, run each with `--no-push` (distinct output files, no git race) and make one commit.

There are **two parallel notebooks**, one per dataset, with identical structure — they differ only in the data sources they load (and the example gene):
- **`wormbhuvanfitter.ipynb`** — the worm dataset (`worm.csv` → `master`, `worm_fourparam_table_excluded_at_or_below_-1.csv` → `fourparam_df`). This is the original notebook (renamed from `newbhuvanfitter.ipynb`; the worm CSVs were renamed from `Supplementary Data 1_csv.csv` / `fourparam_table*`).
- **`cerebellumbhuvanfitter.ipynb`** — the GTEx cerebellum dataset (`cerebellumlog2.csv`, loaded with `.drop(columns=['Description']).set_index('Name').T`, → `master`; `cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv` → `fourparam_df`).

Both are the interactive scratch space (`from bhuvanfitter import BhuvanFitter, gene_peaks`) and do three things. When you change one notebook's shared analysis structure, mirror it in the other:
- **Per-gene inspection** — build a `BhuvanFitter` on one column of `master` (the transposed source CSV) and view its `fit("fourparam")` / `fit("kde")` and the `hist(lines=["fourparam", "kde"])` overlay.
- **Cross-gene distribution plots** — it loads that dataset's excluded fourparam table into `fourparam_df` and defines two helpers for histogramming any column across genes: `select(param)` applies the **single shared filter** (`fit_success == True`, `0 < truncationindex < 1`, and `n_obs >= MIN_OBS`, NaNs dropped) and `plot_param_hist(param, *, color, bins, log)` draws the histogram (returning the data) and prints a filtering funnel — gene counts at each stage (master total → fourparam_df total → fit_success → NaN drop → truncationindex > 0 → < 1 → n_obs >= MIN_OBS). Change the filter in one place (`select`) and both the sumsquarevalue and truncationindex plots follow. **`MIN_OBS` is a notebook-level constant (default 30)** applied as an analysis-time floor on the table's `n_obs` column — the generator's own `MIN_OBS=10` is too low for a 4-param fit to a 40-bin histogram, and filtering here is identical to regenerating the table with a higher floor (a gene's fit is independent of `MIN_OBS` once `n_obs` clears it), so no regeneration is needed to tune it.
- **Per-gene expression summaries (filtered)** — `plot_gene_stat_hist(stat, source='full', *, color, bins)` is a third shared helper (alongside `select` / `plot_param_hist`) that histograms, over the **same filtered genes** (`select('gene')`), a per-gene **mean** or **std** of expression. `source='full'` reduces `master` over **all** samples (incl. the −1 floor); `source='table'` instead plots the post-exclusion `fourparam_df[stat]` column (only `> −1` samples) — genuinely different numbers, so pair the mean/std calls with the **same `source`** if feeding the scatter. The last three cells are now `strain_avgs = plot_gene_stat_hist('mean')`, `strain_sds = plot_gene_stat_hist('std')`, plus the unchanged **mean-vs-SD scatter** (`strain_avgs` x, `strain_sds` y; aligned by gene index). Run them in order (the scatter reuses `strain_avgs` / `strain_sds`).
- **Genes-of-interest truncation ranges (worm only)** — `wormbhuvanfitter.ipynb` has a final section that loads `genes_of_interest.json` (keys `mco_dev` / `mco_behavior` / `lof_dev` / `lof_behavior`, each a gene→transcript-ID map whose IDs match the table's `gene` column), builds **7 categories** (the 4 keys + deduped-union `mco_combined` / `lof_combined` / `all_combined`), and reports per-category `truncationindex` summaries (`n_requested`, `n_in_table`, `n_used`, min/max/mean/median) via `ti_ranges(apply_filter)` in **two cells**: one without the shared filter (any successfully-fit transcript) and one with it. Worm-specific (the transcript IDs are worm) — no cerebellum counterpart.
- **pTriplo vs. fourparam features (cerebellum only)** — `cerebellumbhuvanfitter.ipynb` has a final section ("5. pTriplo vs. the fourparam features") that correlates **pTriplo** (Collins 2022 triplosensitivity score, `dosage_sensitivity_collins2022.tsv.gz`; "high" cutoff ≥ 0.94) against **every numeric fourparam feature**. pTriplo and the truncation index are both *continuous classifier scores* for triplosensitivity (neither is ground truth); this is the descriptive agreement view. The fourparam table is keyed by Ensembl IDs and pTriplo by gene symbols, so the join uses the GTEx source's own two id columns (`Name` = ENSG, `Description` = symbol — version-stripped) as the ENSG→symbol map; no external lookup. Restricted to `fit_success == True` (≈15.5k genes joined), it prints `ptriplo_corr` (Spearman ρ + Pearson r + n per feature, sorted by |ρ|) and a `pTriplo` vs. `feat` scatter (default `feat="truncationindex"`). Cerebellum-specific (pTriplo is human; the worm table's transcript IDs would need ortholog mapping) — no worm counterpart. **Result: `truncationindex` is ~uncorrelated with pTriplo (ρ ≈ −0.01); the expression-level features (`mean`/`min`/`x0`/`max`, ρ ≈ 0.43) carry the signal.**
- **Truncation index conditioned on fit quality (cerebellum only)** — section "6. Truncation index vs dosage scores, conditioned on fit quality" tests whether `truncationindex` only tracks dosage sensitivity for genes the Gaussian actually fits. Fit quality is **R² = 1 − sumsquarevalue/SST** (scale-free), *not* raw `sumsquarevalue` (which is entangled with n_obs/std/A and would bias the correlation). **No refit** — `y0,A,x0,w` and `sumsquarevalue` come from the table; each gene is only re-histogrammed via `BhuvanFitter(vals>−1, gene)` (the constructor caches `hist_counts`; no `curve_fit`) to get `SST = Σ(counts − mean)²`. Builds `cond` (the `select()`-filtered genes inner-joined to **both `pTriplo` and `pHaplo`** + the `r2` column), then (a) a **stratified trend** — Spearman ρ(`truncationindex`, score) within R² deciles, to see if |ρ| rises with fit quality — and (b) `cond_corr`, the full feature-vs-pTriplo/pHaplo correlation table on the well-fit subset (`r2 >= R2_MIN`, default 0.9), plus a pTriplo-vs-truncationindex scatter colored by R². **Result: the hypothesis is NOT supported** — ρ(`truncationindex`, pTriplo/pHaplo) stays ≈ 0 across all R² deciles and on the well-fit subset (≈0.03 / 0.07 at R² ≥ 0.9), while expression-level features keep tracking the scores. Cerebellum-specific.
- **Predicting pTriplo from distribution features (cerebellum only)** — section "7. Predicting pTriplo from distribution features" asks the multivariate question: can a model reconstruct pTriplo from a gene's distribution **shape**? Two feature lists: `LEVEL` (mean/min/max/x0/std/n_obs/A) and scale-free `SHAPE` (skew, kurt, r2, cv=std/mean, w_x0=w/x0, truncationindex, ti_fourparam_sigma_dist). `build_shape_features(df)` adds these with **no refit** (skew/kurt from the post-exclusion values via `scipy.stats`, `r2` from `BhuvanFitter`'s cached histogram, ratios). It runs a **falsifiable EDA gate** (univariate ρ per feature + KMeans clustering of `SHAPE` → per-cluster pTriplo, a "target group" check), a **5-fold CV increment table** (`level`/`shape`/`combined` × Ridge/HistGBM; the deliverable is whether `combined` beats `level`), then fits the combined HistGBM, prints permutation importance, plots out-of-fold predicted-vs-actual, and **applies the model genome-wide** to every `fit_success` gene → `cerebellum_shape_ptriplo_pred.csv` (`gene, symbol, shape_pTriplo, pTriplo`; ~47.8k genes). `TARGET` is a cell-level constant (`"pTriplo"`; set to `"pHaplo"` for the identical pipeline). `n_peaks` was deliberately **not** used (KDE ~0.18 s/gene → infeasible genome-wide, and it added no signal). Requires `scikit-learn`. **Result: shape-only ≈ level-only (ρ ≈ 0.36 vs 0.38, R² ≈ 0.14 vs 0.15) and combined ≈ level — shape and level are redundant, not additive, with a low ceiling (R² ≈ 0.15)**; single-tissue distribution features reconstruct only a fraction of pTriplo. Cerebellum-specific (the trained model could later be applied to the worm table).

## Other data files

`Supplementary Data 1 trunc 20250702.xlsx` (gene-name ↔ identifier mapping) is an input for downstream gene-of-interest analysis. `genes_of_interest.json` (curated worm gene sets: `mco_dev` / `mco_behavior` / `lof_dev` / `lof_behavior`, each gene→transcript-ID) **is consumed by `wormbhuvanfitter.ipynb`'s final section** (truncation-index ranges per category). `dosage_sensitivity_collins2022.tsv.gz` (Collins 2022 gnomAD dosage-sensitivity scores: columns `#gene`, `pHaplo`, `pTriplo`; human gene symbols; ~18.6k genes) **is consumed by `cerebellumbhuvanfitter.ipynb`'s final section** (pTriplo↔fourparam-feature correlation). High-dosage-sensitivity cutoffs are **pTriplo ≥ 0.94** (triplosensitive) and **pHaplo ≥ 0.86** (haploinsufficient) — different thresholds, don't conflate them. `positiveANDnegativeControlGenes.csv` is a list of ~727 human genes found as an extra copy in more than one healthy person (duplication-tolerant / OE-tolerant negative controls). None of these are consumed by the generator scripts.

## The grant proposal (`ACFROG~1.PDF` → `acfrog.md`)

`ACFROG~1.PDF` is the Pierce-lab NIH-style grant proposal ("OE threshold / Censoring Hypothesis" for DS & AD genes) whose preliminary figures this whole toolkit supports. `acfrog.md` is its markdown conversion (text reconstructed from the PDF; Figs 1 & 2 are conceptual diagrams/tables and are only described — but Fig 2A's gene table is transcribed). **`regenerate_acfrog_figures.py`** regenerates the two *data* figures **faithfully from repo data** via `BhuvanFitter` → `acfrog_figure3_worm.png` (Fig 3, worm) and `acfrog_figure4_cerebellum.png` (Fig 4, cerebellum: Hsa21 & triplication-sensitive gene distributions). Run: `python regenerate_acfrog_figures.py`.

- **Fig 3 OE groups come from Figure 2A, not `genes_of_interest.json`.** The script hard-codes `SENS_GENES` (24 OE-sensitive) / `TOL_GENES` (23 OE-tolerant) worm gene names extracted from Fig 2A's mcOE-"Any" shaded cells (the json lacks the OE-tolerant no-phenotype genes, which are the correct "absent" control for panel D). `build_worm_groups()` resolves those names → fourparam-table transcript IDs (`w{n}_{transcript}`) via the `Supplementary Data 1 trunc 20250702.xlsx` `GeneName`/`transcript`/`wwww` columns (hlh-34, sod-5, H39E23.3 aren't in the dataset and drop out). **`EXCLUDE_GENES = {"itsn-1", "adr-2"}`** (two OE-tolerant genes) are dropped from `SENS_GENES`/`TOL_GENES` at the professor's instruction, so they're absent from every worm panel (3A–D and the 3B comparison), filtered and no-filter alike.
- **Honesty caveat (in script + md):** with the proper Fig 2A groups, Fig 3D shows the qualitative pattern as a **tail effect** — OE-tolerant genes are capped (max TI ≈ 0.30) while OE-sensitive spread to 0.60 — but the medians are similar so Mann-Whitney is only p ≈ 0.24 (n.s.), weaker than the proposal's p < 0.004 (likely the older *unbounded* TI metric amplified it). Cerebellum (Fig 4): triplication genes have truncation index ≈ 0, and even as a pTriplo≥0.94-vs-dup-tolerant contrast the signal is flat/reversed (consistent with the section-5 ρ ≈ −0.01). Fig 3D uses **probability (per-bin fraction) normalization** to match the original's `p` axis, not density.
- **`figure3(exclude=...)` toggle.** `figure3(exclude=True)` (default) uses the `> -1` excluded table + floor → `acfrog_figure3_worm.png`; `figure3(exclude=False)` uses the non-excluded `worm_fourparam_table.csv` and includes the −1 floor when re-histogramming → `acfrog_figure3_worm_nofilter.png` (comparison only; the validity filters `fit_success`/`0<TI<1`/`n_obs≥30` still apply). Dropping the exclusion lets more transcripts through and shifts fits left but leaves the qualitative result unchanged. `__main__` writes both.
- **`figure3b_comparison(exclude=...)`**: side-by-side of panel 3B as (left) a box plot over all transcripts (median) vs (right) the **grant's style** — one canonical transcript per gene, mean ± SEM, gene-name labels, Welch t-test. `exclude=True` → `acfrog_figure3b_comparison.png` (filtered: box MWU p ≈ 0.22, grant-style p ≈ 0.11); `exclude=False` → `acfrog_figure3b_comparison_nofilter.png` (no-filter: p ≈ 0.20 / ≈ 0.08). Each figure's suptitle states which. Both exclude itsn-1/adr-2 (dropping those two elevated tolerant genes moved the grant-style p from ≈ 0.29 → ≈ 0.11). `__main__` writes both. Same metric/data/groups; the grant's stronger *p < 0.004* is **not** reproduced. Diagnosis (also in acfrog.md): the box-vs-grant difference is **transcript selection + mean-vs-median statistic**, NOT the `> -1` filter (no-filter barely moves p) and NOT the truncationindex formula (current metric ≈ the grant's own supplementary `trunc` for most genes: cle-1/dnsn-1/dip-2/Y54E10A.11/Y105E8A.1 all within ~0.05). Even the grant's own `trunc` on these 23/21 gene lists only reaches p ≈ 0.03–0.07, so the published p likely used an expanded gene set (+13 independent-study genes).
