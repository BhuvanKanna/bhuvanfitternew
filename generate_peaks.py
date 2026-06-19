#!/usr/bin/env python
"""
generate_peaks.py

Detect the KDE density peaks of every gene in ``worm.csv``,
build a nested dictionary of them, write it to ``peaks.json``, and push that file
to the GitHub repo.

Output structure
----------------
    {
      gene_id: {
        peak_expression_value: {"height": <kde density>, "prominence": <prom>},
        ...
      },
      ...
    }

The number of peaks for a gene is ``len()`` of its inner dict. Genes with no
detectable mode (degenerate data) map to an empty dict.

Note: JSON object keys are strings, so the peak expression values become string
keys in ``peaks.json`` — parse them back to float on load.

Peak detection (``gene_peaks``) and the CSV loader / git push helpers are reused
from the existing modules; this script only orchestrates the per-gene loop.

Usage
-----
    python generate_peaks.py            # all genes, write peaks.json, push
    python generate_peaks.py --no-push  # write peaks.json but skip the git push
    python generate_peaks.py --limit 50 # only the first 50 genes (testing)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from bhuvanfitter import gene_peaks
from generate_fourparam_stats import load_expression, git_push

HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE / "worm.csv"
OUTPUT_JSON = HERE / "peaks.json"


def build_peaks(df, limit: int | None = None) -> dict:
    """Detect peaks for every gene (column) and collect the nested dict."""
    genes = list(df.columns)
    if limit is not None:
        genes = genes[:limit]

    peaks_by_gene: dict = {}
    total = len(genes)
    for i, gene in enumerate(genes, 1):
        peaks_by_gene[gene] = gene_peaks(df[gene].astype(float).values)
        if i % 2000 == 0 or i == total:
            print(f"  peaks {i}/{total} genes")
    return peaks_by_gene


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-push", action="store_true",
                        help="Write peaks.json but do not commit/push to GitHub.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only process the first N genes (for quick testing).")
    args = parser.parse_args()

    print(f"Reading {INPUT_CSV.name} ...")
    df = load_expression(INPUT_CSV)
    print(f"  {df.shape[1]} genes x {df.shape[0]} strains")

    print("Detecting peaks ...")
    peaks_by_gene = build_peaks(df, limit=args.limit)
    counts = [len(p) for p in peaks_by_gene.values()]
    n_with_peaks = sum(1 for c in counts if c > 0)
    print(f"  {n_with_peaks}/{len(peaks_by_gene)} genes have >=1 peak; "
          f"total peaks = {sum(counts)}")

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(peaks_by_gene, f, indent=1)
    print(f"Wrote {OUTPUT_JSON.name} ({len(peaks_by_gene)} genes)")

    if args.no_push:
        print("--no-push set; skipping git push.")
        return

    git_push(
        HERE, OUTPUT_JSON,
        "Update peaks.json (regenerated from worm.csv)",
    )


if __name__ == "__main__":
    main()
