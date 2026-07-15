#!/usr/bin/env python
"""Test the professor's confound hypothesis via Enrichr: are POS/TOL (and our own
truncation-index "ceiling"/"non-ceiling" groups) confounded with tissue-of-origin,
specifically brain vs. non-brain expression -- rather than genuinely reflecting
overexpression tolerance?

Uses Enrichr's REST API (maayanlab.cloud/Enrichr) -- rigorous hypergeometric/Fisher's
exact enrichment testing against curated tissue-expression gene-set libraries, not our
own ad hoc GTEx-median ratio counting (gene_tissue_prominence.csv). Confirmed current
library names for every category the professor named: GTEx_Tissues_V8_2023 (GTEx),
ARCHS4_Tissues, Human_Gene_Atlas (closest available match -- there is no library
literally named "Human_Protein_Atlas" in Enrichr today), Jensen_TISSUES, and
Allen_Brain_Atlas_up (Allen Brain Atlas).

Four gene lists tested:
  POS          - positive_genes_compiled.txt (70 genes)
  TOL          - positiveANDnegativeControlGenes.csv (~648-677 genes)
  ceiling      - cerebellum genes with truncationindex > 0.3 (this repo's existing
                 "clear ceiling" star threshold from the grant figures work)
  non_ceiling  - cerebellum genes with truncationindex == 0 (fully uncapped)

Run: python enrichr_tissue_enrichment.py
"""
import json
import subprocess
import time

import numpy as np
import pandas as pd

from regenerate_grant_figures import load_tol
from train_dosage_classifier import TABLE, symbol_map

BASE = "https://maayanlab.cloud/Enrichr"


def _curl_json(args):
    """Shell out to curl (requests' TLS handshake times out in this environment,
    curl works reliably against the same host)."""
    result = subprocess.run(["curl", "-s", "--max-time", "60"] + args,
                            capture_output=True, text=True, check=True)
    return json.loads(result.stdout)
LIBRARIES = ["GTEx_Tissues_V8_2023", "ARCHS4_Tissues", "Human_Gene_Atlas",
            "Jensen_TISSUES", "Allen_Brain_Atlas_up"]
BRAIN_KEYWORDS = ["brain", "cerebell", "cortex", "hippocamp", "amygdala",
                 "hypothalamus", "striatum", "thalamus", "cerebral", "neuron",
                 "spinal cord", "pons", "medulla", "basal ganglia"]

SUMMARY_OUT = "outputs/tables/enrichr_tissue_summary.csv"
FIG_OUT = "outputs/figures/enrichr_tissue_comparison.png"


def submit_list(genes, description):
    """Writes the gene list to a temp file and has curl read it as the form field's
    content (large lists like TOL/ceiling would risk the OS command-line length
    limit if inlined directly as a -F argument)."""
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("\n".join(genes))
        path = f.name
    try:
        data = _curl_json(["-X", "POST", f"{BASE}/addList",
                           "-F", f"list=<{path}",
                           "-F", f"description={description}"])
    finally:
        import os
        os.unlink(path)
    return data["userListId"]


def get_enrichment(user_list_id, library):
    data = _curl_json([f"{BASE}/enrich?userListId={user_list_id}&backgroundType={library}"])
    rows = data[library]
    cols = ["rank", "term", "pvalue", "zscore", "combined_score", "genes",
           "adj_pvalue", "old_pvalue", "old_adj_pvalue"]
    return pd.DataFrame(rows, columns=cols)


def is_brain_term(term):
    t = term.lower()
    return any(k in t for k in BRAIN_KEYWORDS)


def top_brain_hit(df, top_n=50):
    """Best (lowest adj_pvalue) brain-related term within the top_n ranked terms."""
    sub = df[df["rank"] <= top_n].copy()
    sub = sub[sub["term"].apply(is_brain_term)]
    if len(sub) == 0:
        return None
    return sub.sort_values("adj_pvalue").iloc[0]


def build_ceiling_groups(sample_size=650, seed=0):
    """Cerebellum genes with truncationindex > 0.1 (ceiling, 586 genes genome-wide --
    a meaningfully-sized "has a clear ceiling" group, comparable to TOL's 677) vs a
    random sample of truncationindex == 0 genes (non-ceiling, "fully uncapped").
    TI==0 alone is 27,667 genes genome-wide (the known over-represented uncapped pile
    this repo's own select()/valid() filters normally exclude) -- far too large for one
    Enrichr call and not a fair size-matched comparison, so it's downsampled here."""
    sym = symbol_map()
    df = pd.read_csv(TABLE)
    df = df[(df["fit_success"] == True) & (df["n_obs"] >= 30)].copy()  # noqa: E712
    df["symbol"] = df["gene"].str.split(".").str[0].map(sym)
    df = df.dropna(subset=["symbol"])
    df = df.sort_values("n_obs", ascending=False).drop_duplicates("symbol", keep="first")

    ceiling = sorted(df.loc[df["truncationindex"] > 0.1, "symbol"].unique())
    non_ceiling_all = df.loc[df["truncationindex"] == 0, "symbol"].unique()
    non_ceiling = sorted(pd.Series(non_ceiling_all).sample(
        n=min(sample_size, len(non_ceiling_all)), random_state=seed))
    return ceiling, non_ceiling


def smoke_test():
    print("== Smoke test: known brain genes (GAD1, SYN1, RBFOX3) vs non-brain (PCSK9) ==")
    uid = submit_list(["GAD1", "SYN1", "RBFOX3"], "brain_smoketest")
    time.sleep(1)
    df = get_enrichment(uid, "GTEx_Tissues_V8_2023")
    hit = top_brain_hit(df, top_n=10)
    print(f"Brain genes top brain term (rank<=10): "
         f"{hit['term'] if hit is not None else 'NONE FOUND -- unexpected'} "
         f"(adj_p={hit['adj_pvalue'] if hit is not None else None})")
    assert hit is not None, "Smoke test failed: known brain genes show no brain enrichment"

    uid2 = submit_list(["PCSK9"], "nonbrain_smoketest")
    time.sleep(1)
    df2 = get_enrichment(uid2, "GTEx_Tissues_V8_2023")
    hit2 = top_brain_hit(df2, top_n=5)
    print(f"PCSK9 top brain term in top 5 (expect none/weak): "
         f"{hit2['term'] if hit2 is not None else 'none (as expected)'}")
    print("Smoke test passed.\n")


def run_list(name, genes, results):
    print(f"\n=== {name} (n={len(genes)}) ===")
    uid = submit_list(genes, name)
    time.sleep(1)
    for lib in LIBRARIES:
        try:
            df = get_enrichment(uid, lib)
        except Exception as e:
            print(f"  {lib}: FAILED ({e})")
            continue
        df.to_csv(f"outputs/tables/enrichr_{name}_{lib}.csv", index=False)
        top = df.iloc[0]
        brain_hit = top_brain_hit(df, top_n=50)
        print(f"  {lib:24s} top term: {top['term'][:50]:50s} adj_p={top['adj_pvalue']:.2e}"
             + (f"  | best brain term (rank {int(brain_hit['rank'])}): "
                f"{brain_hit['term'][:40]} adj_p={brain_hit['adj_pvalue']:.2e}"
                if brain_hit is not None else "  | no brain term in top 50"))
        results.append({
            "list": name, "library": lib, "n_genes": len(genes),
            "top_term": top["term"], "top_term_adj_p": top["adj_pvalue"],
            "brain_term": brain_hit["term"] if brain_hit is not None else None,
            "brain_term_rank": int(brain_hit["rank"]) if brain_hit is not None else None,
            "brain_term_adj_p": brain_hit["adj_pvalue"] if brain_hit is not None else np.nan,
            "brain_term_odds_ratio": brain_hit["combined_score"] if brain_hit is not None else np.nan,
        })
        time.sleep(1)


def main():
    smoke_test()

    pos = sorted({l.strip() for l in open("data/positive_genes_compiled.txt") if l.strip()})
    tol = sorted(load_tol())
    ceiling, non_ceiling = build_ceiling_groups()
    print(f"POS={len(pos)}  TOL={len(tol)}  ceiling={len(ceiling)}  non_ceiling={len(non_ceiling)}")

    results = []
    for name, genes in [("POS", pos), ("TOL", tol),
                        ("ceiling", ceiling), ("non_ceiling", non_ceiling)]:
        run_list(name, genes, results)

    summary = pd.DataFrame(results)
    summary.to_csv(SUMMARY_OUT, index=False)
    print(f"\nSaved summary -> {SUMMARY_OUT}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary["neg_log10_adj_p"] = -np.log10(summary["brain_term_adj_p"].clip(lower=1e-300))
    pivot = summary.pivot(index="library", columns="list", values="neg_log10_adj_p")
    pivot = pivot.reindex(columns=["POS", "TOL", "ceiling", "non_ceiling"])

    fig, ax = plt.subplots(figsize=(10, 6))
    pivot.plot(kind="bar", ax=ax, color=["crimson", "steelblue", "darkorange", "gray"])
    ax.set_ylabel("-log10(adjusted p-value), best brain-related term")
    ax.set_title("Enrichr brain-term enrichment: POS/TOL vs ceiling/non-ceiling")
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=1, label="p=0.05")
    ax.legend(title="")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(FIG_OUT, dpi=110)
    print(f"Saved figure -> {FIG_OUT}")

    print("\n== Verdict ==")
    print(summary[["list", "library", "brain_term_rank", "brain_term_adj_p"]]
         .to_string(index=False))


if __name__ == "__main__":
    main()
