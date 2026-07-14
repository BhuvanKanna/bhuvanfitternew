#!/usr/bin/env python
"""Fit each labelled (POS/TOL) gene's 4-parameter Gaussian in ITS OWN most-relevant
tissue, instead of forcing every gene through cerebellum. Every truncation-index
analysis in this repo so far has used cerebellum data uniformly, which is only
biologically relevant for a fraction of genes (TOL is heavily testis/immune-
restricted, established via the Enrichr work) -- analyzing a testis-specific
gene's distribution shape in cerebellum, where it's barely expressed, is unlikely
to reveal a meaningful "ceiling" regardless of its true biology.

Data: the full per-sample GTEx TPM matrix (gtex_full_persample_tpm.gct.gz,
59,033 genes x 19,616 samples across all tissues, 2.15GB -- NOT committed to
git, see .gitignore) plus the sample->tissue attributes file
(gtex_sample_attributes.txt). Each gene's own top tissue comes from the
already-computed gtex_tissue_specificity.csv (via train_tissue_aware_classifier
.load_tissue_features(), same canonical-row dedup used everywhere else in this
repo). Values are log2(TPM+1)-transformed (NOT reusing cerebellumlog2.csv's -1
placeholder convention -- that was an artifact of how that file was originally
pre-processed; log2(TPM+1) is the standard transform and needs no floor value
since TPM=0 maps cleanly to log2(1)=0).

Run: python generate_own_tissue_fourparam.py
"""
import gzip

import numpy as np
import pandas as pd

from bhuvanfitter import BhuvanFitter
from generate_fourparam_stats import COLUMNS, MIN_OBS, _failed_row
from regenerate_grant_figures import load_tol

TPM_MATRIX = "gtex_full_persample_tpm.gct.gz"
SAMPLE_ATTRS = "gtex_sample_attributes.txt"
TISSUE_SPECIFICITY = "gtex_tissue_specificity.csv"
POS_FILE = "positive_genes_compiled.txt"
OUT_TABLE = "gene_own_tissue_fourparam.csv"


def load_top_tissue_by_symbol():
    """symbol -> top_tissue (raw name, not the numeric TISSUE feature block) --
    same canonical-row dedup (highest top_tissue_tpm) used everywhere else."""
    t = pd.read_csv(TISSUE_SPECIFICITY)
    t = t.sort_values("top_tissue_tpm", ascending=False).drop_duplicates("symbol", keep="first")
    return t.set_index("symbol")["top_tissue"]

# Known typo in GTEx's own median-TPM-file column name ('Lymphode' vs 'Lymphoid').
TISSUE_NAME_OVERRIDES = {
    "Small_Intestine_Terminal_Ileum_Lymphode_Aggregate":
        "Small Intestine - Terminal Ileum - Lymphoid Aggregate",
}


def normalize_smtsd(s):
    return s.replace(" - ", "_").replace("(", "").replace(")", "").replace(" ", "_")


def build_bridge(median_tissue_cols, smtsd_values):
    """median-file tissue column name -> raw SMTSD string."""
    norm_smtsd = {normalize_smtsd(s): s for s in smtsd_values}
    bridge = {}
    for col in median_tissue_cols:
        if col in norm_smtsd:
            bridge[col] = norm_smtsd[col]
        elif col in TISSUE_NAME_OVERRIDES:
            bridge[col] = TISSUE_NAME_OVERRIDES[col]
        else:
            raise ValueError(f"Could not bridge tissue column: {col}")
    return bridge


def main():
    top_tissue_by_symbol = load_top_tissue_by_symbol()  # indexed by symbol
    pos = {l.strip() for l in open(POS_FILE) if l.strip()}
    tol = set(load_tol())
    labelled = {s: 1 for s in pos} | {s: 0 for s in tol if s not in pos}

    gene_top_tissue = {s: top_tissue_by_symbol.loc[s] for s in labelled
                       if s in top_tissue_by_symbol.index}
    print(f"Labelled genes: {len(labelled)}  with a tissue-table match: {len(gene_top_tissue)}")

    # Bridge SMTSD (sample attributes) <-> median-file tissue column naming.
    with gzip.open("gtex_median_tpm_by_tissue.gct.gz", "rt") as f:
        f.readline(); f.readline()
        median_cols = f.readline().strip().split("\t")[2:]
    attrs = pd.read_csv(SAMPLE_ATTRS, sep="\t", usecols=["SAMPID", "SMTSD"])
    bridge = build_bridge(median_cols, attrs["SMTSD"].dropna().unique())
    print(f"Tissue bridge covers all {len(bridge)} median-file tissue columns.")

    sampid_to_smtsd = dict(zip(attrs["SAMPID"], attrs["SMTSD"]))
    reverse_bridge = {smtsd: col for col, smtsd in bridge.items()}

    # Read the big matrix's header once to map sample columns -> bridged tissue name.
    with gzip.open(TPM_MATRIX, "rt") as f:
        f.readline(); f.readline()
        header = f.readline().rstrip("\n").split("\t")
    sample_cols = header[2:]
    sample_tissue = [reverse_bridge.get(sampid_to_smtsd.get(s)) for s in sample_cols]

    needed_tissues = set(gene_top_tissue.values())
    tissue_colidx = {t: np.array([i for i, tt in enumerate(sample_tissue) if tt == t])
                     for t in needed_tissues}
    for t in sorted(needed_tissues):
        print(f"  {t}: {len(tissue_colidx[t])} samples available")

    wanted_symbols = set(gene_top_tissue)
    records = []
    found = set()
    with gzip.open(TPM_MATRIX, "rt") as f:
        f.readline(); f.readline(); f.readline()  # header already parsed above
        for line_i, line in enumerate(f, 1):
            parts = line.rstrip("\n").split("\t")
            symbol = parts[1]
            if symbol not in wanted_symbols or symbol in found:
                continue
            found.add(symbol)
            top_tissue = gene_top_tissue[symbol]
            idx = tissue_colidx[top_tissue]
            values = np.array([float(parts[2 + i]) for i in idx], dtype=float)
            log2v = np.log2(values + 1)
            n_obs = int(np.isfinite(log2v).sum())
            log2v = log2v[np.isfinite(log2v)]

            if n_obs < MIN_OBS:
                row = _failed_row(symbol, n_obs)
            else:
                try:
                    bf = BhuvanFitter(log2v, gene_name=symbol)
                    row = bf.fit("fourparam")
                except RuntimeError:
                    row = _failed_row(symbol, n_obs)
            row["tissue"] = top_tissue
            row["n_tissue_samples"] = len(idx)
            records.append(row)

            if len(found) % 100 == 0:
                print(f"  fit {len(found)}/{len(wanted_symbols)} genes "
                     f"(scanned {line_i} matrix rows)")
            if len(found) == len(wanted_symbols):
                break

    missing = wanted_symbols - found
    if missing:
        print(f"Not found in the TPM matrix (symbol mismatch): {sorted(missing)}")

    out = pd.DataFrame.from_records(records, columns=COLUMNS + ["tissue", "n_tissue_samples"])
    out["label"] = out["gene"].map(labelled)
    out.to_csv(OUT_TABLE, index=False)
    print(f"\nWrote {len(out)} genes -> {OUT_TABLE}")
    print(f"fit_success: {out['fit_success'].sum()} / {len(out)}")


if __name__ == "__main__":
    main()
