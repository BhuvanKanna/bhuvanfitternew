#!/usr/bin/env python
"""
generate_fourparam_stats.py

Fit the 4-parameter Gaussian to every gene in ``worm.csv``,
write the per-gene results to ``worm_fourparam_table.csv``, and push that file to the
GitHub repo.

The ``BhuvanFitter`` class is imported from ``bhuvanfitter.py``, the single
source of truth for the fitting logic. The output table has one row per gene
with exactly the columns ``fit("fourparam")`` returns.

This is the **unfiltered** generator: every finite expression value is kept
(no ``<= threshold`` exclusion -- see ``generate_fourparam_stats_excluded.py``
for that variant). It is dataset-agnostic via ``--input`` / ``--id-col`` /
``--drop-col`` (the shared ``load_expression``); the output filename is derived
from the input stem as ``<stem>_fourparam_table.csv`` so datasets never collide
(worm -> ``worm_fourparam_table.csv``, GTEx cerebellum ->
``cerebellumlog2_v8_fourparam_table.csv``).

Usage
-----
    python generate_fourparam_stats.py            # fit all worm genes, write CSV, push
    python generate_fourparam_stats.py --no-push  # write CSV but skip the git push
    python generate_fourparam_stats.py --limit 50 # only the first 50 genes (testing)
    python generate_fourparam_stats.py --input data/cerebellumlog2_v8.csv \
        --id-col Name --drop-col Description --jobs 11 --no-push  # GTEx cerebellum
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from bhuvanfitter import BhuvanFitter

HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE / "data/worm.csv"
OUTPUT_CSV = HERE / "outputs/tables/worm_fourparam_table.csv"


def output_csv_for(input_path: Path) -> Path:
    """Output path whose name is derived from the input dataset stem:
    ``<input stem>_fourparam_table.csv`` (worm -> worm_fourparam_table.csv,
    cerebellumlog2_v8 -> cerebellumlog2_v8_fourparam_table.csv)."""
    return HERE / f"outputs/tables/{input_path.stem}_fourparam_table.csv"

# The exact columns returned by BhuvanFitter.fit("fourparam"), in order.
COLUMNS = [
    "gene", "y0", "A", "x0", "w", "sumsquarevalue",
    "ti_fourparam_sigma_dist", "truncationindex",
    "min", "max", "mean", "std", "skew", "kurt", "right", "maxheight", "rightheight",
    "n_obs", "fit_success",
]
MIN_OBS = 10  # genes with fewer finite observations are not fit (matches the notebook)


def load_expression(csv_path: Path, id_col: str = "strain",
                    drop_cols=()) -> pd.DataFrame:
    """
    Load an expression matrix and orient it as rows = samples, columns = genes
    (the input format BhuvanFitter expects), mirroring the notebook's loader.

    The input file is expected to have one row per gene: ``id_col`` is the
    gene-identifier column (its values become the gene/column labels after the
    transpose) and the remaining columns are samples. ``drop_cols`` lists any
    extra non-sample index columns to discard first (e.g. a ``Description``
    column alongside a ``Name`` id column).

    Parameters
    ----------
    csv_path : Path
        CSV with genes as rows.
    id_col : str
        Column holding the gene identifier. Defaults to ``"strain"`` (the worm
        Supplementary Data 1 layout, whose first column is mislabeled ``strain``
        but actually holds gene names).
    drop_cols : iterable of str
        Additional columns to drop before transposing (non-sample metadata).
    """
    df = pd.read_csv(csv_path)
    if drop_cols:
        df = df.drop(columns=list(drop_cols))
    df = df.set_index(id_col)
    df = df.T  # rows = samples, columns = genes
    return df


def _failed_row(gene: str, n_obs: int) -> dict:
    """A results row for a gene that was skipped or whose fit failed."""
    row = {col: np.nan for col in COLUMNS}
    row["gene"] = gene
    row["n_obs"] = n_obs
    row["fit_success"] = False
    return row


# Set in each worker process by _init_worker so _fit_one can read it without
# re-pickling it per gene.
_WORKER_MAX_NFEV: int = 10_000


def _init_worker(max_nfev: int) -> None:
    global _WORKER_MAX_NFEV
    _WORKER_MAX_NFEV = max_nfev


def _fit_one(item) -> dict:
    """Fit a single gene given ``(gene_name, values_array)``. Module-level and
    picklable so it can run in a multiprocessing pool. Keeps every finite value
    (unfiltered); returns a ``_failed_row`` (never raises) on too-few
    observations or a non-converging fit."""
    gene, values = item
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    n_obs = int(data.size)

    if n_obs < MIN_OBS:
        return _failed_row(gene, n_obs)
    try:
        bf = BhuvanFitter(data, gene_name=gene)
        return bf.fit("fourparam", max_nfev=_WORKER_MAX_NFEV)
    except RuntimeError:
        return _failed_row(gene, n_obs)


def build_table(df: pd.DataFrame, BhuvanFitter, limit: int | None = None,
                max_nfev: int = 10_000, jobs: int = 1) -> pd.DataFrame:
    """Fit every gene (column) and collect one results row per gene.

    ``jobs`` worker processes fit genes in parallel (genes are independent);
    ``jobs=1`` runs in-process (bit-identical, ``imap`` preserves order).
    ``max_nfev`` caps each curve_fit."""
    genes = list(df.columns)
    if limit is not None:
        genes = genes[:limit]
    total = len(genes)
    items = ((gene, df[gene].to_numpy(dtype=float)) for gene in genes)

    records = []
    if jobs <= 1:
        _init_worker(max_nfev)
        for i, gene in enumerate(genes, 1):
            records.append(_fit_one((gene, df[gene].to_numpy(dtype=float))))
            if i % 2000 == 0 or i == total:
                print(f"  fit {i}/{total} genes", flush=True)
    else:
        with mp.Pool(jobs, initializer=_init_worker, initargs=(max_nfev,)) as pool:
            for i, rec in enumerate(pool.imap(_fit_one, items, chunksize=64), 1):
                records.append(rec)
                if i % 2000 == 0 or i == total:
                    print(f"  fit {i}/{total} genes", flush=True)

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
                             "Genes are independent, so this scales near-linearly.")
    parser.add_argument("--max-nfev", type=int, default=10_000, dest="max_nfev",
                        help="curve_fit evaluation cap per gene (default: 10000).")
    args = parser.parse_args()

    input_csv = args.input
    output_csv = output_csv_for(input_csv)

    print(f"Reading {input_csv.name} ...")
    df = load_expression(input_csv, id_col=args.id_col, drop_cols=args.drop_cols)
    print(f"  {df.shape[1]} genes x {df.shape[0]} samples")

    jobs = args.jobs if args.jobs > 0 else (os.cpu_count() or 1)
    print(f"Fitting genes (unfiltered, max_nfev={args.max_nfev}, jobs={jobs}) ...")
    table = build_table(df, BhuvanFitter, limit=args.limit,
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
        f"Update {output_csv.name} (regenerated from {input_csv.name})",
    )


if __name__ == "__main__":
    main()
