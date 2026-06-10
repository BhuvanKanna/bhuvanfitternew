#!/usr/bin/env python
"""
generate_fourparam_stats_excluded.py

Identical to ``generate_fourparam_stats.py`` except each gene's expression
array is first filtered to **drop any value <= -0.75** before fitting. Low
expression values are treated as below-threshold and excluded, so the number of
observations per gene may decrease (and a gene may drop below ``MIN_OBS``).

Fit the 4-parameter Gaussian to every (filtered) gene in
``Supplementary Data 1_csv.csv``, write the per-gene results to
``fourparam_table_excluded.csv``, and push that file to the GitHub repo.

The ``BhuvanFitter`` class is imported from ``bhuvanfitter.py``, the single
source of truth for the fitting logic. The output table has one row per gene
with exactly the columns ``fit("fourparam")`` returns.

Usage
-----
    python generate_fourparam_stats_excluded.py            # fit all genes, write CSV, push
    python generate_fourparam_stats_excluded.py --no-push  # write CSV but skip the git push
    python generate_fourparam_stats_excluded.py --limit 50 # only the first 50 genes (testing)
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
OUTPUT_CSV = HERE / "fourparam_table_excluded.csv"

# Expression values at or below this threshold are excluded before fitting.
EXCLUDE_AT_OR_BELOW = -0.75


def build_table(df: pd.DataFrame, BhuvanFitter, limit: int | None = None) -> pd.DataFrame:
    """Fit every gene (column) after dropping values <= EXCLUDE_AT_OR_BELOW."""
    genes = list(df.columns)
    if limit is not None:
        genes = genes[:limit]

    records = []
    total = len(genes)
    for i, gene in enumerate(genes, 1):
        data = df[gene].astype(float).values
        data = data[np.isfinite(data)]
        data = data[data > EXCLUDE_AT_OR_BELOW]  # drop low/below-threshold expression
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
    args = parser.parse_args()

    print(f"Reading {INPUT_CSV.name} ...")
    df = load_expression(INPUT_CSV)
    print(f"  {df.shape[1]} genes x {df.shape[0]} strains")

    print(f"Fitting genes (excluding values <= {EXCLUDE_AT_OR_BELOW}) ...")
    table = build_table(df, BhuvanFitter, limit=args.limit)
    n_ok = int(table["fit_success"].sum())
    print(f"  {n_ok}/{len(table)} genes fit successfully")

    table.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {OUTPUT_CSV.name} ({len(table)} rows, {len(COLUMNS)} columns)")

    if args.no_push:
        print("--no-push set; skipping git push.")
        return

    git_push(
        HERE, OUTPUT_CSV,
        "Update fourparam_table_excluded.csv (values <= -0.75 excluded before fit)",
    )


if __name__ == "__main__":
    main()
