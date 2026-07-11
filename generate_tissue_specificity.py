#!/usr/bin/env python
"""Genome-wide GTEx tissue-specificity features, to test whether OE-tolerant
genes are less prominently expressed in brain (relative to other tissues) than
OE-sensitive genes.

Source: GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz -- the official
GTEx median-TPM-per-tissue summary matrix (59,033 genes x 68 tissue/sub-tissue
columns, incl. 13 distinct brain sub-regions), downloaded once from
https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/ and committed
to the repo (~8.8MB compressed). Values are linear TPM (unlike cerebellumlog2.csv,
which is already log2) -- log2(TPM+1)-transformed here before any differencing.

Per gene, computes:
  brain_cerebellum          log2 median TPM in Brain_Cerebellum specifically
                            (continuity with this repo's existing focus tissue)
  brain_mean / brain_max    mean/max log2 median TPM across all 13 Brain_* columns
  nonbrain_mean / nonbrain_max  same, across the ~55 non-brain tissues
  log2fc_brain_vs_nonbrain  brain_mean - nonbrain_mean (already log2 scale)
  tau                       standard 0-1 tissue-specificity index (computed on
                            linear TPM across all 68 tissues; tau=1 -> expressed
                            in essentially one tissue, tau=0 -> uniform)
  top_tissue / top_tissue_tpm   name and linear-TPM value of the single highest-
                            median-TPM tissue

Output: gtex_tissue_specificity.csv (ensg, symbol, + the above), genome-wide.

Run: python generate_tissue_specificity.py
"""
import gzip
import shutil
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

URL = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
      "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
GCT_GZ = Path("gtex_median_tpm_by_tissue.gct.gz")
OUT = "gtex_tissue_specificity.csv"

BRAIN_COLS = [
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra",
]


def ensure_downloaded():
    if GCT_GZ.exists():
        return
    print(f"Downloading {URL} -> {GCT_GZ}")
    urllib.request.urlretrieve(URL, GCT_GZ)


def load_gct():
    with gzip.open(GCT_GZ, "rt") as f:
        f.readline()  # "#1.2"
        f.readline()  # "<n_genes>\t<n_cols>"
        df = pd.read_csv(f, sep="\t")
    return df


def tau(row_linear):
    """Standard tissue-specificity index (Yanai et al. 2005): 0 = uniform,
    1 = expressed in a single tissue. Computed on linear TPM."""
    m = row_linear.max()
    if m <= 0:
        return np.nan
    xhat = row_linear / m
    return (1 - xhat).sum() / (len(row_linear) - 1)


def main():
    ensure_downloaded()
    df = load_gct()
    tissue_cols = [c for c in df.columns if c not in ("Name", "Description")]
    nonbrain_cols = [c for c in tissue_cols if c not in BRAIN_COLS]
    print(f"Loaded {len(df)} genes x {len(tissue_cols)} tissues "
         f"({len(BRAIN_COLS)} brain, {len(nonbrain_cols)} non-brain)")

    linear = df[tissue_cols].astype(float)
    log2v = np.log2(linear + 1)

    out = pd.DataFrame({
        "ensg": df["Name"].str.split(".").str[0],
        "symbol": df["Description"],
        "brain_cerebellum": log2v["Brain_Cerebellum"],
        "brain_mean": log2v[BRAIN_COLS].mean(axis=1),
        "brain_max": log2v[BRAIN_COLS].max(axis=1),
        "nonbrain_mean": log2v[nonbrain_cols].mean(axis=1),
        "nonbrain_max": log2v[nonbrain_cols].max(axis=1),
    })
    out["log2fc_brain_vs_nonbrain"] = out["brain_mean"] - out["nonbrain_mean"]
    out["tau"] = linear.apply(tau, axis=1)
    top_idx = linear.values.argmax(axis=1)
    out["top_tissue"] = np.array(tissue_cols)[top_idx]
    out["top_tissue_tpm"] = linear.values[np.arange(len(linear)), top_idx]

    out.to_csv(OUT, index=False)
    print(f"Wrote {len(out)} genes -> {OUT}")
    print(out[["brain_cerebellum", "brain_mean", "nonbrain_mean",
              "log2fc_brain_vs_nonbrain", "tau"]].describe().round(3))


if __name__ == "__main__":
    main()
