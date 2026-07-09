# Grant figures 3B, 3D, 10A driven by the positive-gene list

**Date:** 2026-07-09
**Reference document:** `grant.pdf` (the current NIH proposal; supersedes `acfrog.md` for this work)

## Goal

Regenerate grant **Figure 3B**, **Figure 3D**, and **Figure 10A** using the
user-supplied human dosage-sensitive "positive" gene list (`positive_genes.txt`)
as an OE-sensitive group, **shown alongside** the grant's own OE-sensitive set,
against an OE-tolerant control, in **both** the worm and human (cerebellum)
datasets (10A is human-only, per the grant).

## What the grant figures are (from grant.pdf)

- **Fig 3B** (worm): truncation index (ratio of the truncated right side of the
  expression distribution to the peak) compared between genes whose
  overexpression *caused* an OE phenotype vs. those that did not — "genes that
  caused OE phenotypes had significantly higher indices." A group comparison.
- **Fig 3D** (worm): the Fig-3B data as **normalized histograms** (per-bin
  fraction, `p` axis) for phenotype-present vs. phenotype-absent vs. **all
  transcripts**.
- **Fig 10A** (human, GTEx cerebellum): **cumulative (CDF)** curves of the
  truncation index. All cerebellar genes (63% have trunc=0); segmental-
  duplication / pangenome-tolerant genes shift **left** (OE-tolerant control);
  predicted OE-sensitive genes (pTriplo / OE-sensitive CNV / dominant-GF) shift
  **right**. Groups compared with **KS tests**. The key readout is the fraction
  uncapped (trunc=0): ~0.70 tolerant vs ~0.45 sensitive.

## Gene sets

Each set is resolved in **both** human-symbol space and worm-ortholog space.

| Set | Source | Role |
|-----|--------|------|
| **POS** | `positive_genes.txt` (44 human symbols) | OE-sensitive (user's list) |
| **GRANT** | grant Fig 2A mcOE-phenotype genes | OE-sensitive (grant's list) |
| **TOL** | `positiveANDnegativeControlGenes.csv` (~726 duplication-tolerant symbols) | OE-tolerant control |
| **ALL** | every gene in the dataset | background reference |

- **GRANT** in worm space = the existing `SENS_GENES` worm names (Fig 2A mcOE
  "Any" = present, minus `itsn-1`/`adr-2`). GRANT in human space = the Fig-2A
  human partners of those worm names (e.g. `sod-1`→SOD1, `rcan-1`→RCAN1,
  `D1037.1`→SON, `cle-1`→COL18A1, …). Merged paralog rows (KCNJ6/KCNJ15,
  RRP1/RRP1B, USP25/USP28) contribute all listed human symbols.
- **POS** and **TOL** are already human symbols → used directly in human space;
  mapped to worm orthologs for worm space.
- **TOL parsing:** read column 0 of `positiveANDnegativeControlGenes.csv`, drop
  the `geneDUPtol` header token, keep entries matching a gene-symbol pattern
  (uppercase alnum), dedupe. (The other columns hold free-text citation
  fragments and are ignored.)

## Ortholog bridge (worm versions only)

Map POS + TOL human symbols → *C. elegans* orthologs → worm transcript IDs:

1. Fetch human↔worm orthologs from Alliance/DIOPT (per-gene API, proven to work;
   ~770 genes: 44 POS + 726 TOL). Cache the result to **`human_worm_orthologs.tsv`**
   (`human_symbol, wbgene_id, worm_gene, worm_symbol_source`) so the run is
   reproducible and offline on reruns. Refetch only genes not already cached.
2. worm gene name → worm transcript ID via the existing `build_worm_groups`
   machinery (`Supplementary Data 1 trunc 20250702.xlsx` `GeneName`/`transcript`/
   `wwww` columns → `w{n}_{transcript}`).

GRANT worm IDs come from `build_worm_groups(SENS_GENES)` directly (no ortholog
lookup needed). Genes with no ortholog / not in the worm table simply drop out;
the script reports coverage counts (n_requested → n_mapped → n_in_table → n_valid).

## Truncation-index data & validity filter

- Worm: `worm_fourparam_table_excluded_at_or_below_-1.csv`, `truncationindex`.
- Human: `cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv`,
  `truncationindex`; join the table's Ensembl `gene` (ENSG) to human symbols via
  the GTEx source's `Name`↔`Description` map (version-stripped), exactly as
  `cerebellumbhuvanfitter.ipynb` does.
- **Validity filter (deliberately different from the notebooks):**
  `fit_success == True & truncationindex.notna() & 0 <= truncationindex <= 1 &
  n_obs >= 30`. **TI = 0 (uncapped) genes are KEPT** — they carry the whole
  tolerant-vs-sensitive signal and the grant explicitly counts them (the notebook
  `select`/`valid` uses `0 < TI < 1` and must NOT be reused here).
- A gene set maps to a set of transcript rows; per-figure each set's TI values are
  the valid rows for those transcripts (worm: all isoforms; human: one row/ENSG).

## Figures produced

New standalone script **`regenerate_grant_figures.py`** (leaves
`regenerate_acfrog_figures.py` untouched). Reuses `build_worm_groups`,
`curve_xy`, `BhuvanFitter`, and the constants from the existing script (import or
copy the small pieces; no behavior change to the original).

Each figure overlays POS, GRANT, TOL, and ALL:

1. `grant_figure3b_worm.png` / `grant_figure3b_human.png` — box/strip plot of
   `truncationindex` per group, with pairwise **Mann-Whitney U** p-values (POS vs
   TOL, GRANT vs TOL, POS vs GRANT) and per-group n in labels.
2. `grant_figure3d_worm.png` / `grant_figure3d_human.png` — per-bin-fraction
   normalized histograms (step style, matching the grant's `p` axis) of POS,
   GRANT, TOL overlaid on ALL.
3. `grant_figure10a_human.png` — cumulative (CDF) curves of `truncationindex`
   for POS, GRANT, TOL, ALL, with **KS-test** p-values vs ALL and the
   trunc=0 uncapped fraction annotated per group. Cerebellum-only (grant's 10A).

`__main__` writes all five PNGs and prints the coverage + stats summary.

## Non-goals (YAGNI)

- No worm 10A (grant's 10A is human-only).
- No refit of any gene — read `truncationindex` from the committed tables.
- No changes to `regenerate_acfrog_figures.py`, the generators, or the notebooks.
- No `--no-filter`/exclude toggle variants (the committed excluded tables only).

## Risks / notes

- **Ortholog fetch cost:** first run makes ~770 API calls (minutes). Mitigated by
  caching to `human_worm_orthologs.tsv`; commit the cache so reruns are instant.
- **Worm coverage of TOL/POS may be modest** (many human dosage genes are not
  1:1-conserved). The coverage report makes this explicit; small n is expected
  and acceptable for the worm panels.
- **Expected result** (consistent with the repo's prior findings): the human
  `truncationindex` separates OE-tolerant from OE-sensitive only weakly
  (ρ≈−0.01 vs pTriplo), so POS may not shift strongly from TOL. The figures are
  faithful reproductions of the grant's *method*; whether they reproduce the
  grant's *effect size* is an empirical outcome to report, not to force.
