#!/usr/bin/env python
"""
generate_fourparam_stats_excluded.py

Identical to ``generate_fourparam_stats.py`` except each gene's expression
array is first filtered to **drop any value <= EXCLUDE_AT_OR_BELOW** before
fitting. Low expression values are treated as below-threshold and excluded, so
the number of observations per gene may decrease (and a gene may drop below
``MIN_OBS``).

Fit the 4-parameter Gaussian to every (filtered) gene in ``worm.csv``, write the
per-gene results to a CSV whose name **encodes the input dataset and the exclusion
threshold** -- ``<input stem>_fourparam_table_excluded_at_or_below_<threshold>.csv``
(e.g. ``worm_fourparam_table_excluded_at_or_below_-1.csv``) -- and push that file to
the GitHub repo. Because the output name is derived from the input and
``EXCLUDE_AT_OR_BELOW``, changing the threshold (or ``--input``) automatically writes
a separate, self-labeled spreadsheet instead of overwriting another one.

The ``BhuvanFitter`` class is imported from ``bhuvanfitter.py``, the single
source of truth for the fitting logic. The output table has one row per gene
with exactly the columns ``fit("fourparam")`` returns.

Usage
-----
    python generate_fourparam_stats_excluded.py                  # fit all genes (threshold -1), write CSV, push
    python generate_fourparam_stats_excluded.py --no-push        # write CSV but skip the git push
    python generate_fourparam_stats_excluded.py --limit 50       # only the first 50 genes (testing)
    python generate_fourparam_stats_excluded.py --threshold -0.75  # use a different exclusion threshold
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import BhuvanFitter
from generate_fourparam_stats import COLUMNS, MIN_OBS, load_expression, git_push, _failed_row

HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE / "worm.csv"

# Expression values at or below this threshold are excluded before fitting.
# This is the default; the --threshold CLI flag overrides it per run.
EXCLUDE_AT_OR_BELOW = -1


def output_csv_for(threshold, input_path: Path = INPUT_CSV) -> Path:
    """Output path whose name encodes the **input dataset** and the exclusion
    threshold, so each setting writes its own spreadsheet and datasets never
    collide: ``<input stem>_fourparam_table_excluded_at_or_below_<threshold>.csv``
    (e.g. worm -> worm_fourparam_table_excluded_at_or_below_-1.csv, cerebellumlog2
    -> cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv). ``:g`` keeps
    -1.0 rendered as ``-1`` so int and float thresholds map to the same name."""
    return HERE / f"{input_path.stem}_fourparam_table_excluded_at_or_below_{threshold:g}.csv"


# Default output (threshold = EXCLUDE_AT_OR_BELOW); --threshold overrides it.
OUTPUT_CSV = output_csv_for(EXCLUDE_AT_OR_BELOW)


# Set in each worker process by _init_worker so _fit_one can read them without
# re-pickling them per gene.
_WORKER_THRESHOLD: float = EXCLUDE_AT_OR_BELOW
_WORKER_MAX_NFEV: int = 10_000


def _init_worker(exclude_at_or_below: float, max_nfev: int) -> None:
    global _WORKER_THRESHOLD, _WORKER_MAX_NFEV
    _WORKER_THRESHOLD = exclude_at_or_below
    _WORKER_MAX_NFEV = max_nfev


def _fit_one(item) -> dict:
    """Fit a single gene given ``(gene_name, values_array)``. Module-level and
    picklable so it can run in a multiprocessing pool. Drops values
    <= the worker's threshold first; returns a ``_failed_row`` (never raises) on
    too-few observations or a non-converging fit."""
    gene, values = item
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    data = data[data > _WORKER_THRESHOLD]  # drop low/below-threshold expression
    n_obs = int(data.size)

    if n_obs < MIN_OBS:
        return _failed_row(gene, n_obs)
    try:
        bf = BhuvanFitter(data, gene_name=gene)
        return bf.fit("fourparam", max_nfev=_WORKER_MAX_NFEV)
    except RuntimeError:
        return _failed_row(gene, n_obs)


def build_table(df: pd.DataFrame, BhuvanFitter, limit: int | None = None,
                exclude_at_or_below: float = EXCLUDE_AT_OR_BELOW,
                max_nfev: int = 10_000, jobs: int = 1) -> pd.DataFrame:
    """Fit every gene (column) after dropping values <= exclude_at_or_below.

    ``jobs`` worker processes fit genes in parallel (genes are independent);
    ``jobs=1`` runs in-process. ``max_nfev`` caps each curve_fit."""
    genes = list(df.columns)
    if limit is not None:
        genes = genes[:limit]
    total = len(genes)
    # Stream (gene, values) pairs so we never materialise all arrays at once.
    items = ((gene, df[gene].to_numpy(dtype=float)) for gene in genes)

    records = []
    if jobs <= 1:
        _init_worker(exclude_at_or_below, max_nfev)
        for i, gene in enumerate(genes, 1):
            records.append(_fit_one((gene, df[gene].to_numpy(dtype=float))))
            if i % 2000 == 0 or i == total:
                print(f"  fit {i}/{total} genes", flush=True)
    else:
        with mp.Pool(jobs, initializer=_init_worker,
                     initargs=(exclude_at_or_below, max_nfev)) as pool:
            for i, rec in enumerate(pool.imap(_fit_one, items, chunksize=64), 1):
                records.append(rec)
                if i % 2000 == 0 or i == total:
                    print(f"  fit {i}/{total} genes", flush=True)

    return pd.DataFrame.from_records(records, columns=COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-push", action="store_true",
                        help="Write the CSV but do not commit/push to GitHub.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only fit the first N genes (for quick testing).")
    parser.add_argument("--threshold", type=float, default=EXCLUDE_AT_OR_BELOW,
                        help="Exclude expression values <= this threshold before "
                             "fitting (default: %(default)s). The output filename "
                             "encodes it: fourparam_table_excluded_at_or_below_<threshold>.csv.")
    parser.add_argument("--input", type=Path, default=INPUT_CSV,
                        help="Input CSV (genes as rows). Defaults to the worm "
                             "Supplementary Data 1 file. A non-default input gets "
                             "its stem prefixed onto the output filename.")
    parser.add_argument("--id-col", type=str, default="strain",
                        help="Gene-identifier column in the input (default: strain "
                             "for the worm file; use e.g. Name for the GTEx data).")
    parser.add_argument("--drop-col", action="append", default=[], dest="drop_cols",
                        metavar="COL",
                        help="Extra non-sample column to drop before transposing "
                             "(repeatable, e.g. --drop-col Description).")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Worker processes to fit genes in parallel (default: 1). "
                             "Genes are independent, so this scales near-linearly; "
                             "e.g. --jobs 11 on a 12-core machine.")
    parser.add_argument("--max-nfev", type=int, default=10_000, dest="max_nfev",
                        help="curve_fit evaluation cap per gene (default: 10000). "
                             "Lowering it (e.g. 2000) mainly speeds up non-converging "
                             "genes, which otherwise burn the whole budget.")
    args = parser.parse_args()

    threshold = args.threshold
    input_csv = args.input
    output_csv = output_csv_for(threshold, input_csv)

    print(f"Reading {input_csv.name} ...")
    df = load_expression(input_csv, id_col=args.id_col, drop_cols=args.drop_cols)
    print(f"  {df.shape[1]} genes x {df.shape[0]} samples")

    jobs = args.jobs if args.jobs > 0 else (os.cpu_count() or 1)
    print(f"Fitting genes (excluding values <= {threshold:g}, "
          f"max_nfev={args.max_nfev}, jobs={jobs}) ...")
    table = build_table(df, BhuvanFitter, limit=args.limit,
                        exclude_at_or_below=threshold,
                        max_nfev=args.max_nfev, jobs=jobs)
    n_ok = int(table["fit_success"].sum())
    print(f"  {n_ok}/{len(table)} genes fit successfully")

    table.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv.name} ({len(table)} rows, {len(COLUMNS)} columns)")

    if args.no_push:
        print("--no-push set; skipping git push.")
        return

    git_push(
        HERE, output_csv,
        f"Update {output_csv.name} (values <= {threshold:g} excluded before fit)",
    )


if __name__ == "__main__":
    main()
