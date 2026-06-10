#!/usr/bin/env python
"""
generate_fourparam_stats.py

Fit the 4-parameter Gaussian to every gene in ``Supplementary Data 1_csv.csv``,
write the per-gene results to ``fourparam_table.csv``, and push that file to the
GitHub repo.

The ``BhuvanFitter`` class is imported from ``bhuvanfitter.py``, the single
source of truth for the fitting logic. The output table has one row per gene
with exactly the columns ``fit("fourparam")`` returns.

Usage
-----
    python generate_fourparam_stats.py            # fit all genes, write CSV, push
    python generate_fourparam_stats.py --no-push  # write CSV but skip the git push
    python generate_fourparam_stats.py --limit 50 # only the first 50 genes (testing)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import BhuvanFitter

HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE / "Supplementary Data 1_csv.csv"
OUTPUT_CSV = HERE / "fourparam_table.csv"

# The exact columns returned by BhuvanFitter.fit("fourparam"), in order.
COLUMNS = [
    "gene", "y0", "A", "x0", "w", "sumsquarevalue",
    "ti_fourparam_sigma_dist", "truncationindex",
    "min", "max", "right", "maxheight", "rightheight",
    "n_obs", "fit_success",
]
MIN_OBS = 10  # genes with fewer finite observations are not fit (matches the notebook)


def load_expression(csv_path: Path) -> pd.DataFrame:
    """
    Load the expression matrix and orient it as rows = strains, columns = genes
    (the input format BhuvanFitter expects), mirroring the notebook's loader.
    """
    df = pd.read_csv(csv_path)
    df = df.set_index("strain")
    df = df.T  # rows = strains (isolates), columns = genes
    return df


def _failed_row(gene: str, n_obs: int) -> dict:
    """A results row for a gene that was skipped or whose fit failed."""
    row = {col: np.nan for col in COLUMNS}
    row["gene"] = gene
    row["n_obs"] = n_obs
    row["fit_success"] = False
    return row


def build_table(df: pd.DataFrame, BhuvanFitter, limit: int | None = None) -> pd.DataFrame:
    """Fit every gene (column) and collect one results row per gene."""
    genes = list(df.columns)
    if limit is not None:
        genes = genes[:limit]

    records = []
    total = len(genes)
    for i, gene in enumerate(genes, 1):
        data = df[gene].astype(float).values
        data = data[np.isfinite(data)]
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


def git_push(repo_dir: Path, file_path: Path, message: str) -> None:
    """Stage, commit, and push the output CSV to origin. No-op if unchanged."""
    def run(*args: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(repo_dir), *args],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"`git {' '.join(args)}` failed:\n{result.stderr.strip()}")
        return result.stdout

    run("add", file_path.name)
    if not run("status", "--porcelain", file_path.name).strip():
        print(f"No changes to {file_path.name} — nothing to commit or push.")
        return
    run("commit", "-m", message)
    run("push", "origin", "HEAD")
    print(f"Pushed {file_path.name} to origin.")


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

    print("Fitting genes ...")
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
        "Update fourparam_table.csv (regenerated from Supplementary Data 1_csv.csv)",
    )


if __name__ == "__main__":
    main()
