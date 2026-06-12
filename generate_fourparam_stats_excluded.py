#!/usr/bin/env python
"""
generate_fourparam_stats_excluded.py

Identical to ``generate_fourparam_stats.py`` except each gene's expression
array is first filtered to **drop any value <= EXCLUDE_AT_OR_BELOW** before
fitting. Low expression values are treated as below-threshold and excluded, so
the number of observations per gene may decrease (and a gene may drop below
``MIN_OBS``).

Fit the 4-parameter Gaussian to every (filtered) gene in
``Supplementary Data 1_csv.csv``, write the per-gene results to a CSV whose name
**encodes the exclusion threshold** -- ``fourparam_table_excluded_at_or_below_<threshold>.csv``
(e.g. ``fourparam_table_excluded_at_or_below_-1.csv``) -- and push that file to
the GitHub repo. Because the output name is derived from ``EXCLUDE_AT_OR_BELOW``,
changing the threshold automatically writes a separate, self-labeled spreadsheet
instead of overwriting a differently-thresholded one.

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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import BhuvanFitter
from generate_fourparam_stats import COLUMNS, MIN_OBS, load_expression, git_push, _failed_row

HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE / "Supplementary Data 1_csv.csv"

# Expression values at or below this threshold are excluded before fitting.
# This is the default; the --threshold CLI flag overrides it per run.
EXCLUDE_AT_OR_BELOW = -1


def output_csv_for(threshold) -> Path:
    """Output path whose name encodes the exclusion threshold, so each setting
    writes its own spreadsheet (e.g. -1 -> fourparam_table_excluded_at_or_below_-1.csv,
    -0.75 -> fourparam_table_excluded_at_or_below_-0.75.csv). ``:g`` keeps -1.0
    rendered as ``-1`` so int and float thresholds map to the same name."""
    return HERE / f"fourparam_table_excluded_at_or_below_{threshold:g}.csv"


# Default output (threshold = EXCLUDE_AT_OR_BELOW); --threshold overrides it.
OUTPUT_CSV = output_csv_for(EXCLUDE_AT_OR_BELOW)


def build_table(df: pd.DataFrame, BhuvanFitter, limit: int | None = None,
                exclude_at_or_below: float = EXCLUDE_AT_OR_BELOW) -> pd.DataFrame:
    """Fit every gene (column) after dropping values <= exclude_at_or_below."""
    genes = list(df.columns)
    if limit is not None:
        genes = genes[:limit]

    records = []
    total = len(genes)
    for i, gene in enumerate(genes, 1):
        data = df[gene].astype(float).values
        data = data[np.isfinite(data)]
        data = data[data > exclude_at_or_below]  # drop low/below-threshold expression
        n_obs = int(data.size)

        if n_obs < MIN_OBS:
            records.append(_failed_row(gene, n_obs))
        else:
            try:
                bf = BhuvanFitter(data, gene_name=gene)
                records.append(bf.fit("fourparam"))
            except RuntimeError as exc:
                print(f"  warn: fit failed for '{gene}': {exc}", file=sys.stderr)
                records.append(_failed_row(gene, n_obs))

        if i % 2000 == 0 or i == total:
            print(f"  fit {i}/{total} genes")

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
    args = parser.parse_args()

    threshold = args.threshold
    output_csv = output_csv_for(threshold)

    print(f"Reading {INPUT_CSV.name} ...")
    df = load_expression(INPUT_CSV)
    print(f"  {df.shape[1]} genes x {df.shape[0]} strains")

    print(f"Fitting genes (excluding values <= {threshold:g}) ...")
    table = build_table(df, BhuvanFitter, limit=args.limit, exclude_at_or_below=threshold)
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
