"""Regenerate grant figures 3B, 3D, 10A driven by the positive-gene list.

Reference: grant.pdf (Figs 3B, 3D, 10A). See
docs/superpowers/specs/2026-07-09-grant-figures-positives-design.md.

Overlays POS (positive_genes.txt) + GRANT (Fig-2A mcOE set) vs TOL
(positiveANDnegativeControlGenes.csv, duplication-tolerant control) vs ALL, in
both the worm and human (GTEx cerebellum) datasets. 10A is human-only (as in the
grant). Truncation-index values are read from the committed excluded fourparam
tables (no refit). Human<->worm orthologs are fetched once from HGNC + Alliance
and cached to human_worm_orthologs.tsv.

Validity filter: fit_success & 0 < truncationindex < 1 & n_obs >= 30. The open
interval OMITS the truncationindex == 0 (uncapped) and == 1 (fully-capped) piles,
which are over-represented and likely partly artefactual, at the user's request.
(This matches the notebooks' / regenerate_acfrog_figures.py `0 < TI < 1` filter;
an earlier version of this script kept TI == 0 to count the uncapped fraction.)
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

POS_TXT = "data/positive_genes.txt"
TOL_CSV = "data/positiveANDnegativeControlGenes.csv"
ORTH_CACHE = "data/human_worm_orthologs.tsv"

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


# ---------------------------------------------------------------------------
# human -> worm orthologs (HGNC symbol -> HGNC id -> Alliance -> C. elegans)
# ---------------------------------------------------------------------------
def _http_json(url, accept="application/json"):
    req = urllib.request.Request(
        url, headers={"User-Agent": "curl/8", "Accept": accept})
    return json.load(urllib.request.urlopen(req, timeout=30))


def _alliance_worm_orthologs(human_symbol):
    """Human gene symbol -> list of C. elegans ortholog gene symbols (Alliance/
    DIOPT). Resolves the symbol to an HGNC id first (Alliance keys on HGNC)."""
    docs = _http_json(
        f"https://rest.genenames.org/fetch/symbol/{human_symbol}"
    )["response"]["docs"]
    if not docs or "hgnc_id" not in docs[0]:
        return []
    hgnc = docs[0]["hgnc_id"]
    d = _http_json(
        f"https://www.alliancegenome.org/api/gene/{hgnc}/orthologs?stringencyFilter=all")
    out = []
    for r in d.get("results", []):
        o = r.get("geneToGeneOrthologyGenerated", {}).get("objectGene", {})
        if o.get("taxon", {}).get("name") == "Caenorhabditis elegans":
            s = o.get("geneSymbol", {}).get("displayText")
            if s:
                out.append(s)
    return out


def human_to_worm(symbols):
    """Map human symbols -> list of C. elegans ortholog gene symbols, backed by
    the human_worm_orthologs.tsv cache (columns human_symbol, worm_symbol; a
    human with no ortholog gets one row with an empty worm_symbol). Only symbols
    absent from the cache are fetched; the cache is rewritten afterwards."""
    symbols = set(symbols)
    cache = {}
    try:
        cdf = pd.read_csv(ORTH_CACHE, sep="\t").fillna("")
        for h, sub in cdf.groupby("human_symbol"):
            cache[h] = [w for w in sub["worm_symbol"].astype(str) if w]
    except FileNotFoundError:
        pass
    missing = sorted(s for s in symbols if s not in cache)
    for i, s in enumerate(missing):
        try:
            cache[s] = _alliance_worm_orthologs(s)
        except Exception as e:
            print(f"    ortholog fetch failed for {s}: {e}")
            cache[s] = []
        time.sleep(0.1)
        if i % 25 == 0:
            print(f"  orthologs {i}/{len(missing)}")
    rows = []
    for h, ws in sorted(cache.items()):
        if ws:
            rows += [{"human_symbol": h, "worm_symbol": w} for w in ws]
        else:
            rows.append({"human_symbol": h, "worm_symbol": ""})
    pd.DataFrame(rows).to_csv(ORTH_CACHE, sep="\t", index=False)
    return {s: cache.get(s, []) for s in symbols}


# ---------------------------------------------------------------------------
# truncation-index data + group -> TI resolvers
# ---------------------------------------------------------------------------
def valid(table, filtered=True):
    """Validity filter. `filtered=True` uses the OPEN interval 0 < TI < 1
    (drops the over-represented TI==0 uncapped / TI==1 fully-capped piles; the
    notebooks' / acfrog `select`/`valid`, the user's standard process).
    `filtered=False` uses the CLOSED interval 0 <= TI <= 1 (keeps the 0/1 piles;
    faithful to the grant's own Fig-10A method, which counts the uncapped
    fraction). Both also require fit_success and n_obs >= MIN_OBS."""
    ti = table["truncationindex"]
    lo = (ti > 0) if filtered else (ti >= 0)
    hi = (ti < 1) if filtered else (ti <= 1)
    m = table["fit_success"] & ti.notna() & lo & hi & (table["n_obs"] >= MIN_OBS)
    return table[m]


_worm_tab = None
_cereb = None  # (raw cerebellum table indexed by ENSG, symbol -> [ENSG] map)


def _worm_table():
    global _worm_tab
    if _worm_tab is None:
        _worm_tab = pd.read_csv(WORM_TABLE).set_index("gene")
    return _worm_tab


def _cereb_table():
    """Cache the RAW cerebellum table (indexed by ENSG) + the symbol->[ENSG] map;
    the validity filter is applied per call so both filtered/unfiltered modes
    share one read."""
    global _cereb
    if _cereb is None:
        t = pd.read_csv(CEREB_TABLE)
        t["ensg"] = t["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        src = pd.read_csv(CEREB_CSV, usecols=["Name", "Description"])
        src["ensg"] = src["Name"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        sym2ensg = {}
        for _, r in src.iterrows():
            sym2ensg.setdefault(str(r["Description"]), []).append(r["ensg"])
        _cereb = (t.set_index("ensg"), sym2ensg)
    return _cereb


def build_worm_groups_names(names):
    """Resolve arbitrary worm gene names -> fourparam transcript IDs present in
    the worm table (mirrors regenerate_acfrog_figures.build_worm_groups, but for
    any name list). Returns (ids, unmapped_names)."""
    tab = _worm_table()
    xl = pd.read_excel(WORM_MAP_XLSX).dropna(subset=["transcript"])
    xl["tid"] = ("w" + xl["wwww"].astype(int).astype(str) + "_"
                 + xl["transcript"].astype(str))
    xl["seqname"] = xl["transcript"].str.replace(r"\.\d+$", "", regex=True)
    tabset = set(tab.index)
    by_name, by_seq = {}, {}
    for _, r in xl.iterrows():
        if r["tid"] not in tabset:
            continue
        by_name.setdefault(str(r["GeneName"]), []).append(r["tid"])
        by_seq.setdefault(str(r["seqname"]), []).append(r["tid"])
    ids, unmapped = [], []
    for n in names:
        got = (by_name.get(n) or by_seq.get(n)
               or by_seq.get(re.sub(r"\.\d+$", "", str(n))))
        if got:
            ids += got
        else:
            unmapped.append(n)
    return list(dict.fromkeys(ids)), unmapped


def worm_all_ti(filtered=True):
    return valid(_worm_table(), filtered)["truncationindex"].to_numpy()


def human_all_ti(filtered=True):
    raw, _ = _cereb_table()
    return valid(raw, filtered)["truncationindex"].to_numpy()


def worm_ti(worm_names, filtered=True):
    """TI values for the transcripts of the given worm gene names."""
    v = valid(_worm_table(), filtered)
    ids, _ = build_worm_groups_names(worm_names)
    ids = [i for i in ids if i in v.index]
    return v.loc[ids, "truncationindex"].to_numpy() if ids else np.array([])


def human_ti(symbols, filtered=True):
    """TI values for the given human symbols (one row per matched ENSG)."""
    raw, sym2ensg = _cereb_table()
    v = valid(raw, filtered)
    ensgs = [e for s in symbols for e in sym2ensg.get(s, []) if e in v.index]
    ensgs = list(dict.fromkeys(ensgs))
    return v.loc[ensgs, "truncationindex"].to_numpy() if ensgs else np.array([])


# ---------------------------------------------------------------------------
# group assembler
# ---------------------------------------------------------------------------
def groups_for(dataset, filtered=True):
    """Assemble {POS, GRANT, TOL, ALL} -> TI arrays for `dataset` in
    {"worm","human"} at the given filter mode. Prints a per-group coverage line."""
    pos, tol = load_pos(), load_tol()
    if dataset == "human":
        g = {"POS": human_ti(pos, filtered), "GRANT": human_ti(GRANT_HUMAN, filtered),
             "TOL": human_ti(tol, filtered), "ALL": human_all_ti(filtered)}
    elif dataset == "worm":
        h2w = human_to_worm(pos | tol)
        pos_w = [w for s in pos for w in h2w.get(s, [])]
        tol_w = [w for s in tol for w in h2w.get(s, [])]
        g = {"POS": worm_ti(pos_w, filtered), "GRANT": worm_ti(GRANT_WORM, filtered),
             "TOL": worm_ti(tol_w, filtered), "ALL": worm_all_ti(filtered)}
    else:
        raise ValueError(dataset)
    tag = "filtered 0<TI<1" if filtered else "unfiltered 0<=TI<=1"
    for k, v in g.items():
        if len(v):
            print(f"  [{dataset}/{tag}] {k}: n={len(v)}  median={np.median(v):.3f}  "
                  f"frac0={np.mean(v == 0):.2f}")
        else:
            print(f"  [{dataset}/{tag}] {k}: n=0")
    return g


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
COLORS = {"POS": "#d1495b", "GRANT": "#e3a008", "TOL": "#2e86ab", "ALL": "0.55"}
LABELS = {"POS": "POS (positives)", "GRANT": "GRANT (Fig 2A)",
          "TOL": "TOL (dup-tolerant)", "ALL": "All genes"}


def _mwu(a, b):
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    return mannwhitneyu(a, b, alternative="two-sided").pvalue


def _out(name, filtered):
    """Prefix the filtered (0<TI<1) renderings with 'filtered_', and place the
    result under outputs/figures/."""
    return "outputs/figures/" + ("filtered_" if filtered else "") + name


def _tag(filtered):
    return "0<TI<1" if filtered else "0<=TI<=1 (keeps 0/1)"


def figure3b(dataset, filtered=True):
    """Fig 3B: truncation index, POS vs TOL only -- jittered dots with a
    black mean +/- SEM marker per group, Mann-Whitney U."""
    g = groups_for(dataset, filtered)
    order = ["POS", "TOL"]
    fig, ax = plt.subplots(figsize=(5, 5))
    rng = np.random.default_rng(0)
    for i, k in enumerate(order, 1):
        y = g[k]
        if not len(y):
            continue
        ax.scatter(rng.normal(i, 0.06, len(y)), y, s=10, color=COLORS[k],
                   alpha=0.5, linewidths=0, zorder=1)
        m = float(np.mean(y))
        sem = float(np.std(y, ddof=1) / np.sqrt(len(y))) if len(y) > 1 else 0.0
        ax.errorbar(i, m, yerr=sem, fmt="o", color="black", ecolor="black",
                    capsize=6, markersize=7, elinewidth=1.5, zorder=3,
                    label="mean ± SEM" if i == 1 else None)
    ax.set_xlim(0.5, 2.5)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"{k}\n(n={len(g[k])})" for k in order])
    ax.set_ylabel("truncation index")
    p = _mwu(g["POS"], g["TOL"])
    mp, mt = float(np.mean(g["POS"])), float(np.mean(g["TOL"]))
    ax.set_title(f"Fig 3B ({dataset}, {_tag(filtered)}): truncation index, POS vs TOL\n"
                 f"mean±SEM  POS={mp:.3f}  TOL={mt:.3f} | MWU p={p:.3g}",
                 fontsize=9)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    out = _out(f"grant_figure3b_{dataset}.png", filtered)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)
    return out


def _prob_hist(ax, vals, color, label, bins=20):
    if len(vals) == 0:
        return
    counts, edges = np.histogram(vals, bins=bins, range=(0, 1))
    frac = counts / counts.sum()
    centers = (edges[:-1] + edges[1:]) / 2
    ax.step(centers, frac, where="mid", color=color, lw=2,
            label=f"{label} (n={len(vals)})")


def figure3d(dataset, filtered=True):
    """Fig 3D: per-bin-fraction normalized truncation-index histograms."""
    g = groups_for(dataset, filtered)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for k in ["ALL", "TOL", "GRANT", "POS"]:
        _prob_hist(ax, g[k], COLORS[k], LABELS[k])
    ax.set_xlabel("truncation index")
    ax.set_ylabel("fraction of genes (p)")
    ax.set_title(f"Fig 3D ({dataset}, {_tag(filtered)}): "
                 f"normalized truncation-index histograms", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = _out(f"grant_figure3d_{dataset}.png", filtered)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)
    return out


def _cdf(ax, vals, color, label, filtered):
    if len(vals) == 0:
        return
    x = np.sort(vals)
    y = np.arange(1, len(x) + 1) / len(x)
    stat = (f"median={np.median(vals):.3f}" if filtered
            else f"frac0={np.mean(vals == 0):.2f}")
    ax.plot(x, y, color=color, lw=2, label=f"{label} (n={len(vals)}, {stat})")


def figure10a(filtered=True):
    """Fig 10A (human cerebellum): cumulative truncation index + KS vs ALL."""
    g = groups_for("human", filtered)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    _cdf(ax, g["ALL"], "black", LABELS["ALL"], filtered)
    _cdf(ax, g["TOL"], COLORS["TOL"], LABELS["TOL"], filtered)
    _cdf(ax, g["GRANT"], COLORS["GRANT"], LABELS["GRANT"], filtered)
    _cdf(ax, g["POS"], COLORS["POS"], LABELS["POS"], filtered)

    def ks(a):
        return float("nan") if len(a) < 3 else ks_2samp(a, g["ALL"]).pvalue

    ax.set_xlabel("truncation index")
    ax.set_ylabel("cumulative fraction")
    ax.set_title(f"Fig 10A (human cerebellum, {_tag(filtered)}): cumulative truncation index\n"
                 f"KS vs ALL  POS p={ks(g['POS']):.3g} | GRANT p={ks(g['GRANT']):.3g} | "
                 f"TOL p={ks(g['TOL']):.3g}", fontsize=9)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    out = _out("grant_figure10a_human.png", filtered)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)
    return out


if __name__ == "__main__":
    # Two renderings per figure: unfiltered (0<=TI<=1, keeps the 0/1 piles;
    # grant-faithful) and filtered (0<TI<1, 'filtered_' prefix; the notebooks'
    # standard). 10 PNGs total.
    for filt in (False, True):
        for ds in ("human", "worm"):
            figure3b(ds, filt)
            figure3d(ds, filt)
        figure10a(filt)
    print("done: 10 figures (5 unfiltered + 5 filtered)")
