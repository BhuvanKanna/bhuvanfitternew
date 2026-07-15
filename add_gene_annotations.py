#!/usr/bin/env python
"""
add_gene_annotations.py

One-off enrichment: add ``wormbasegeneid`` and ``genename`` columns (right
after ``gene``) to every fourparam table, worm and cerebellum alike, so each
table is self-describing without a separate join.

Worm tables (``gene`` = ``w{n}_{transcript}``):
    ``wormbasegeneid`` / ``genename`` come directly from the
    ``Supplementary Data 1 trunc 20250702.xlsx`` map (``wwww`` + ``transcript``
    reconstruct the ``w{n}_{transcript}`` id, same join key
    ``regenerate_acfrog_figures.build_worm_groups`` already uses;
    ``WormBaseGeneID`` / ``GeneName`` columns supply the values). Full coverage
    for every transcript present in the xlsx.

Cerebellum table (``gene`` = versioned ENSG id):
    ``genename`` comes straight from ``cerebellumlog2.csv``'s own
    ``Name``/``Description`` columns (full coverage — every gene has a GTEx
    symbol). ``wormbasegeneid`` requires crossing species: ENSG -> human
    symbol (``Description``) -> worm ortholog symbol (via the existing
    ``human_worm_orthologs.tsv`` HGNC/Alliance-DIOPT cache) -> WormBaseGeneID
    (via the same xlsx map, keyed by ``GeneName``). That cache only covers the
    ~721 human symbols looked up for the POS/TOL classifier work, so
    ``wormbasegeneid`` is populated only for genes in that set — NaN
    elsewhere, not a genome-wide ortholog fetch.

Run: python add_gene_annotations.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
WORM_MAP_XLSX = HERE / "Supplementary Data 1 trunc 20250702.xlsx"
CEREB_SOURCE = HERE / "cerebellumlog2.csv"
ORTHOLOG_TSV = HERE / "human_worm_orthologs.tsv"

WORM_TABLES = [
    HERE / "worm_fourparam_table.csv",
    HERE / "worm_fourparam_table_excluded_at_or_below_-1.csv",
    HERE / "worm_fourparam_table_excluded_at_or_below_-0.75.csv",
]
CEREB_TABLE = HERE / "cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv"


def _insert_after_gene(df: pd.DataFrame, wormbasegeneid, genename) -> pd.DataFrame:
    df = df.drop(columns=["wormbasegeneid", "genename"], errors="ignore")
    df.insert(1, "genename", genename)
    df.insert(1, "wormbasegeneid", wormbasegeneid)
    return df


def worm_id_map() -> pd.DataFrame:
    xl = pd.read_excel(WORM_MAP_XLSX).dropna(subset=["transcript"])
    xl["tid"] = "w" + xl["wwww"].astype(int).astype(str) + "_" + xl["transcript"].astype(str)
    return xl.drop_duplicates("tid").set_index("tid")[["WormBaseGeneID", "GeneName"]]


def annotate_worm_tables() -> None:
    id_map = worm_id_map()
    for path in WORM_TABLES:
        if not path.exists():
            print(f"skip (missing): {path.name}")
            continue
        df = pd.read_csv(path)
        joined = df[["gene"]].join(id_map, on="gene")
        df = _insert_after_gene(df, joined["WormBaseGeneID"], joined["GeneName"])
        df.to_csv(path, index=False)
        n = joined["WormBaseGeneID"].notna().sum()
        print(f"{path.name}: {n}/{len(df)} rows annotated")


def annotate_cerebellum_table() -> None:
    src = pd.read_csv(CEREB_SOURCE, usecols=["Name", "Description"])
    ensg_to_symbol = src.drop_duplicates("Name").set_index("Name")["Description"]

    orth = pd.read_csv(ORTHOLOG_TSV, sep="\t").dropna(subset=["worm_symbol"])
    orth = orth[orth["worm_symbol"].str.strip() != ""]
    symbol_to_worm_symbol = orth.drop_duplicates("human_symbol").set_index("human_symbol")["worm_symbol"]

    xl = pd.read_excel(WORM_MAP_XLSX).dropna(subset=["GeneName"])
    worm_symbol_to_wbgene = xl.drop_duplicates("GeneName").set_index("GeneName")["WormBaseGeneID"]

    df = pd.read_csv(CEREB_TABLE)
    genename = df["gene"].map(ensg_to_symbol)
    wormbasegeneid = genename.map(symbol_to_worm_symbol).map(worm_symbol_to_wbgene)
    df = _insert_after_gene(df, wormbasegeneid, genename)
    df.to_csv(CEREB_TABLE, index=False)
    print(f"{CEREB_TABLE.name}: genename {genename.notna().sum()}/{len(df)}, "
          f"wormbasegeneid {wormbasegeneid.notna().sum()}/{len(df)} rows annotated")


if __name__ == "__main__":
    annotate_worm_tables()
    annotate_cerebellum_table()
