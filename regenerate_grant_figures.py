"""Regenerate grant figures 3B, 3D, 10A driven by the positive-gene list.

Reference: grant.pdf (Figs 3B, 3D, 10A). See
docs/superpowers/specs/2026-07-09-grant-figures-positives-design.md.

Overlays POS (positive_genes.txt) + GRANT (Fig-2A mcOE set) vs TOL
(positiveANDnegativeControlGenes.csv, duplication-tolerant control) vs ALL, in
both the worm and human (GTEx cerebellum) datasets. 10A is human-only (as in the
grant). Truncation-index values are read from the committed excluded fourparam
tables (no refit). Human<->worm orthologs are fetched once from HGNC + Alliance
and cached to human_worm_orthologs.tsv.

Unlike the notebooks / regenerate_acfrog_figures.py, the validity filter here
KEEPS truncationindex == 0 (uncapped) genes -- they carry the whole
tolerant-vs-sensitive signal and the grant explicitly counts them.
"""
import json
import re
import time
import urllib.request

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ks_2samp

from regenerate_acfrog_figures import (
    MIN_OBS, EXCLUDE_AT_OR_BELOW, WORM_TABLE, CEREB_TABLE, CEREB_CSV,
    WORM_MAP_XLSX, SENS_GENES,
)

POS_TXT = "positive_genes.txt"
TOL_CSV = "positiveANDnegativeControlGenes.csv"
ORTH_CACHE = "human_worm_orthologs.tsv"

# GRANT (Fig-2A mcOE-phenotype) worm names -> human symbol(s), transcribed from
# grant.pdf Figure 2A. Merged paralog rows contribute all listed human symbols.
GRANT_WORM_TO_HUMAN = {
    "chaf-2": ["CHAF1B"], "cle-1": ["COL18A1"], "dip-2": ["DIP2A"],
    "dnsn-1": ["DONSON"], "eva-1": ["EVA1C"], "pat-3": ["ITGB2"],
    "irk-2": ["KCNJ6", "KCNJ15"], "Y54E10A.11": ["LTN1"], "mrps-6": ["MRPS6"],
    "ncam-1": ["NCAM2"], "F43G9.12": ["PAXBP1"], "pdxk-1": ["PDXK"],
    "pfk-1.1": ["PFKL"], "rcan-1": ["RCAN1"], "rrp-1": ["RRP1", "RRP1B"],
    "nrd-1": ["SCAF4"], "Y105E8A.1": ["SH3BGR"], "hlh-34": ["SIM2"],
    "sod-1": ["SOD1"], "D1037.1": ["SON"], "unc-26": ["SYNJ1"],
    "Y74C10AL.2": ["TMEM50B"], "trpp-10": ["TRAPPC10"], "K02C4.3": ["USP25", "USP28"],
}
GRANT_WORM = list(SENS_GENES)  # already itsn-1/adr-2 filtered by the import
GRANT_HUMAN = {h for w in GRANT_WORM for h in GRANT_WORM_TO_HUMAN.get(w, [])}


def load_pos():
    """The user's positive (OE-sensitive) human symbols."""
    return {l.strip() for l in open(POS_TXT) if l.strip()}


def load_tol():
    """Duplication-tolerant (OE-tolerant control) human symbols from column 0 of
    positiveANDnegativeControlGenes.csv (drop the header token + citation text)."""
    df = pd.read_csv(TOL_CSV, header=None, usecols=[0], names=["sym"])
    syms = df["sym"].dropna().astype(str).str.strip()
    ok = syms[syms.str.fullmatch(r"[A-Z0-9][A-Za-z0-9orf\-\.]*")]
    return {s for s in ok if s != "geneDUPtol"}
