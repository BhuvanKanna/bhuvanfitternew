#!/usr/bin/env python
"""Direct answer to "which tissue is each POS/TOL gene most prominent in":
filters gtex_tissue_specificity.csv to the compiled positive genes and the
existing duplication-tolerant control list, for a standalone lookup table
independent of any ML model.

Run: python generate_gene_tissue_prominence.py
"""
import pandas as pd

from regenerate_grant_figures import load_tol

POS_FILE = "positive_genes_compiled.txt"
OUT = "gene_tissue_prominence.csv"


def main():
    tissue = pd.read_csv("gtex_tissue_specificity.csv")
    pos = {l.strip() for l in open(POS_FILE) if l.strip()}
    tol = load_tol()

    rows = []
    for group, symbols in [("POS", pos), ("TOL", tol)]:
        for s in symbols:
            sub = tissue[tissue["symbol"] == s]
            if len(sub) == 0:
                continue
            # one row per symbol: highest top_tissue_tpm among matches (canonical)
            r = sub.sort_values("top_tissue_tpm", ascending=False).iloc[0]
            rows.append({
                "symbol": s, "group": group, "top_tissue": r["top_tissue"],
                "top_tissue_tpm": r["top_tissue_tpm"], "brain_mean": r["brain_mean"],
                "nonbrain_mean": r["nonbrain_mean"],
                "log2fc_brain_vs_nonbrain": r["log2fc_brain_vs_nonbrain"],
                "tau": r["tau"],
            })

    out = pd.DataFrame(rows).sort_values(["group", "log2fc_brain_vs_nonbrain"],
                                         ascending=[True, False])
    out.to_csv(OUT, index=False)
    print(f"Wrote {len(out)} genes -> {OUT} "
         f"(POS={sum(out.group=='POS')}, TOL={sum(out.group=='TOL')})")

    for group in ["POS", "TOL"]:
        g = out[out["group"] == group]
        n_brain_top = g["top_tissue"].str.startswith("Brain").sum()
        print(f"{group}: {n_brain_top}/{len(g)} ({n_brain_top/len(g):.1%}) have a "
             f"Brain_* tissue as their single most-expressed tissue; "
             f"median log2fc_brain_vs_nonbrain={g['log2fc_brain_vs_nonbrain'].median():.3f}, "
             f"median tau={g['tau'].median():.3f}")


if __name__ == "__main__":
    main()
