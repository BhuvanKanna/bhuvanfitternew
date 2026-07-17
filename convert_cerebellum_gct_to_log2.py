#!/usr/bin/env python
"""
convert_cerebellum_gct_to_log2.py

Convert a raw GTEx cerebellum ``gene_tpm`` GCT matrix into this repo's
``cerebellumlog2``-style CSV.

The transform is **``log2(TPM + 1) - 1``**, applied elementwise to every
per-sample TPM value. This is the exact transform the original
``data/cerebellumlog2.csv`` was built with (verified: integer TPMs 0,1,2,3,...
map to -1, 0, 0.585, 1, ... -- so TPM=0 lands on the ``-1`` floor with no
clamping, and every value is >= -1 by construction). The downstream "excluded
at or below -1" fourparam table therefore drops exactly the zero-TPM
(undetected-expression) samples.

Input GCT layout (tab-separated, first two lines are ``#1.3`` and the dims):
``id  Name  Description  <sample columns...>``. The ``id`` column is dropped;
``Name`` (Ensembl id) and ``Description`` (gene symbol) are kept as the two
leading columns, matching ``data/cerebellumlog2.csv``.

Usage
-----
    python convert_cerebellum_gct_to_log2.py \
        --input gene_tpm_brain_cerebellum_v8.gct \
        --output data/cerebellumlog2_v8.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def convert(input_gct: Path, output_csv: Path) -> None:
    print(f"Reading {input_gct.name} ...")
    # Skip the two GCT header lines (`#1.3` and the `<nrow> <ncol> 0 0` dims).
    df = pd.read_csv(input_gct, sep="\t", skiprows=2)
    print(f"  {df.shape[0]} genes x {df.shape[1]} columns (incl. id/Name/Description)")

    if "id" in df.columns:
        df = df.drop(columns=["id"])

    meta = ["Name", "Description"]
    sample_cols = [c for c in df.columns if c not in meta]
    print(f"  transforming {len(sample_cols)} sample columns with log2(TPM+1)-1 ...")

    vals = df[sample_cols].to_numpy(dtype=float)
    vals = np.log2(vals + 1.0) - 1.0
    df[sample_cols] = vals

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Wrote {output_csv} ({df.shape[0]} rows, {df.shape[1]} columns)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path,
                        default=HERE / "gene_tpm_brain_cerebellum_v8.gct",
                        help="Raw GTEx gene_tpm GCT (default: the v8 cerebellum file).")
    parser.add_argument("--output", type=Path,
                        default=HERE / "data/cerebellumlog2_v8.csv",
                        help="Output log2 CSV (default: data/cerebellumlog2_v8.csv).")
    args = parser.parse_args()
    convert(args.input, args.output)


if __name__ == "__main__":
    main()
