#!/usr/bin/env python
"""Rebuild genes_of_interest.json from grant_genes.csv (the full Fig-2A table).

The previous genes_of_interest.json only captured 4 of grant_genes.csv's 6
phenotype columns (mco_dev, mco_behavior, lof_dev, lof_behavior -- missing both
"Any" columns) and had no human-gene names, keyed only by worm gene name (which
collides for irk-2/rrp-1/K02C4.3, each shared by two human paralogs).

This writes a JSON **list** of 49 row-records (one per grant_genes.csv row,
avoiding the worm-gene-name key collisions), each with the human gene, worm
gene, all 6 boolean phenotype flags, and the worm gene's fourparam transcript
IDs (resolved the same way the original JSON's IDs were: via
regenerate_grant_figures.build_worm_groups_names, the xlsx gene<->transcript
mapping).

Run: python update_genes_of_interest_json.py
"""
import csv
import json
import re

import pandas as pd

from regenerate_acfrog_figures import WORM_MAP_XLSX
from regenerate_grant_figures import _worm_table

SRC = "grant_genes.csv"
OUT = "genes_of_interest.json"


def as_bool(cell):
    return cell.strip().lower() == "x"


def build_worm_name_index():
    """(name -> [transcript_ids]) for every worm gene, reading the xlsx ONCE
    (build_worm_groups_names re-reads it per call, far too slow for 48 lookups)."""
    tab = _worm_table()
    xl = pd.read_excel(WORM_MAP_XLSX).dropna(subset=["transcript"])
    xl["tid"] = ("w" + xl["wwww"].astype(int).astype(str) + "_" + xl["transcript"].astype(str))
    xl["seqname"] = xl["transcript"].str.replace(r"\.\d+$", "", regex=True)
    tabset = set(tab.index)
    by_name, by_seq = {}, {}
    for _, r in xl.iterrows():
        if r["tid"] not in tabset:
            continue
        by_name.setdefault(str(r["GeneName"]), []).append(r["tid"])
        by_seq.setdefault(str(r["seqname"]), []).append(r["tid"])

    def lookup(n):
        got = by_name.get(n) or by_seq.get(n) or by_seq.get(re.sub(r"\.\d+$", "", str(n)))
        return list(dict.fromkeys(got)) if got else []
    return lookup


def main():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))

    lookup = build_worm_name_index()
    worm_names = sorted({r["C. elegans"] for r in rows})
    tx_by_worm = {n: lookup(n) for n in worm_names}

    records = []
    for r in rows:
        worm = r["C. elegans"]
        records.append({
            "row": int(r["#"]),
            "human_gene": r["Human Gene"],
            "worm_gene": worm,
            "mcOE_any": as_bool(r["mcOE_Any"]),
            "mcOE_dev": as_bool(r["mcOE_Dev"]),
            "mcOE_behavior": as_bool(r["mcOE_Behavior"]),
            "lof_any": as_bool(r["LOF_Any"]),
            "lof_dev": as_bool(r["LOF_Dev"]),
            "lof_behavior": as_bool(r["LOF_Behavior"]),
            "transcript_ids": tx_by_worm[worm],
        })

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)

    n_mcoe_any = sum(r["mcOE_any"] for r in records)
    n_lof_any = sum(r["lof_any"] for r in records)
    n_with_tx = sum(bool(r["transcript_ids"]) for r in records)
    print(f"Wrote {len(records)} records -> {OUT}")
    print(f"  mcOE_any=True: {n_mcoe_any}   lof_any=True: {n_lof_any}   "
         f"records with >=1 transcript id: {n_with_tx}")


if __name__ == "__main__":
    main()
