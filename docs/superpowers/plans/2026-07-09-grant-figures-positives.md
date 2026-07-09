# Grant Figures 3B/3D/10A (positives) Implementation Plan

> **For agentic workers:** Steps use checkbox (`- [ ]`) syntax. This repo has **no
> pytest/test suite** (per CLAUDE.md); "verification" means running the module or a
> throwaway `python -c` snippet and inspecting printed coverage/stats and that each
> PNG file is written and non-trivial in size. Commit after each task.

**Goal:** Produce five grant figures — `grant_figure3b_worm/human.png`,
`grant_figure3d_worm/human.png`, `grant_figure10a_human.png` — overlaying the user's
positive-gene list (POS) with the grant's Fig-2A OE-sensitive set (GRANT) against a
duplication-tolerant control (TOL) and all genes (ALL), in worm and human data.

**Architecture:** One new standalone script `regenerate_grant_figures.py`. It reuses
`build_worm_groups`, `curve_xy`, constants, and `BhuvanFitter` from the existing
`regenerate_acfrog_figures.py` (imported, not modified). Gene sets resolve to
`truncationindex` values from the committed excluded fourparam tables. Human↔worm
orthologs are fetched once from Alliance/DIOPT and cached to `human_worm_orthologs.tsv`.

**Tech Stack:** Python, numpy, pandas, scipy.stats (mannwhitneyu, ks_2samp),
matplotlib, urllib (Alliance REST), openpyxl (xlsx map).

## Global Constraints

- Truncation-index tables (read-only, no refit): worm `worm_fourparam_table_excluded_at_or_below_-1.csv`; human `cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv`.
- Validity filter (ALL figures): `fit_success == True & truncationindex.notna() & 0 <= truncationindex <= 1 & n_obs >= 30`. **Keep TI==0 genes.**
- MIN_OBS = 30; EXCLUDE_AT_OR_BELOW = -1.0.
- Do NOT modify `regenerate_acfrog_figures.py`, the generators, or the notebooks.
- Human gene keying: table `gene` = ENSG (version-stripped); map ENSG→symbol via the GTEx source `cerebellumlog2.csv` `Name`↔`Description` columns (version-stripped), as `cerebellumbhuvanfitter.ipynb` does.
- GRANT human symbols are the Fig-2A partners of `SENS_GENES` (see Task 1 map).
- Commit + push each task per CLAUDE.md.

---

### Task 1: Gene-set loaders (POS, GRANT-human, TOL)

**Files:**
- Create: `regenerate_grant_figures.py` (header, imports, constants, loaders)

**Interfaces:**
- Produces: `load_pos() -> set[str]`, `load_tol() -> set[str]`, `GRANT_HUMAN: set[str]`, `GRANT_WORM: list[str]` (= `SENS_GENES` imported).

- [ ] **Step 1:** Create the file with imports and reuse from the existing script:
```python
"""Regenerate grant figures 3B, 3D, 10A driven by the positive-gene list.

Reference: grant.pdf (Figs 3B, 3D, 10A). See
docs/superpowers/specs/2026-07-09-grant-figures-positives-design.md.
Overlays POS (positive_genes.txt) + GRANT (Fig-2A mcOE set) vs TOL
(positiveANDnegativeControlGenes.csv) vs ALL, in worm and human (cerebellum).
"""
import json, re, time, urllib.request
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ks_2samp

from regenerate_acfrog_figures import (
    MIN_OBS, EXCLUDE_AT_OR_BELOW, WORM_TABLE, CEREB_TABLE, CEREB_CSV,
    SENS_GENES, build_worm_groups,
)

POS_TXT = "positive_genes.txt"
TOL_CSV = "positiveANDnegativeControlGenes.csv"
ORTH_CACHE = "human_worm_orthologs.tsv"

# GRANT (Fig-2A mcOE-phenotype) worm names -> human symbol(s), from grant.pdf Fig 2A.
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
    return {l.strip() for l in open(POS_TXT) if l.strip()}


def load_tol():
    df = pd.read_csv(TOL_CSV, header=None, usecols=[0], names=["sym"])
    syms = (df["sym"].dropna().astype(str).str.strip())
    ok = syms[syms.str.fullmatch(r"[A-Z0-9][A-Z0-9orf\-\.]*")]
    ok = ok[ok != "geneDUPtol"]
    return set(ok)
```

- [ ] **Step 2:** Verify the loaders:
```bash
python -c "import regenerate_grant_figures as g; print('POS',len(g.load_pos())); print('TOL',len(g.load_tol())); print('GRANT_HUMAN',len(g.GRANT_HUMAN)); print(sorted(g.GRANT_HUMAN)[:5]); print('GRANT_WORM',len(g.GRANT_WORM))"
```
Expected: `POS 44`, `TOL` ~700+, `GRANT_HUMAN` ~26, `GRANT_WORM` 22.

- [ ] **Step 3:** Commit:
```bash
git add regenerate_grant_figures.py && git commit -m "grant figs: gene-set loaders (POS/GRANT/TOL)"
```

---

### Task 2: Human↔worm ortholog fetch + cache

**Files:**
- Modify: `regenerate_grant_figures.py` (add `human_to_worm()`), create `human_worm_orthologs.tsv`

**Interfaces:**
- Consumes: `load_pos`, `load_tol`.
- Produces: `human_to_worm(symbols: set[str]) -> dict[str, list[str]]` (human symbol → list of worm gene symbols), backed by the `human_worm_orthologs.tsv` cache with columns `human_symbol, worm_symbol` (one row per pair; humans with no worm ortholog get one row with empty `worm_symbol`).

- [ ] **Step 1:** Add the fetch+cache function:
```python
def _alliance_worm_orthologs(human_symbol):
    """Return list of C. elegans ortholog gene symbols for a human symbol."""
    url = (f"https://www.alliancegenome.org/api/gene/HGNC:{human_symbol}/orthologs"
           f"?stringencyFilter=all")
    # HGNC id unknown; Alliance also resolves by symbol via search API:
    sr = json.load(urllib.request.urlopen(urllib.request.Request(
        f"https://www.alliancegenome.org/api/search?q={human_symbol}&category=gene&limit=5",
        headers={"User-Agent": "curl/8"}), timeout=30))
    hid = None
    for r in sr.get("results", []):
        if r.get("symbol") == human_symbol and r.get("species") == "Homo sapiens":
            hid = r.get("id"); break
    if not hid:
        return []
    d = json.load(urllib.request.urlopen(urllib.request.Request(
        f"https://www.alliancegenome.org/api/gene/{hid}/orthologs?stringencyFilter=all",
        headers={"User-Agent": "curl/8"}), timeout=30))
    out = []
    for r in d.get("results", []):
        o = r.get("geneToGeneOrthologyGenerated", {}).get("objectGene", {})
        if o.get("taxon", {}).get("name") == "Caenorhabditis elegans":
            s = o.get("geneSymbol", {}).get("displayText")
            if s:
                out.append(s)
    return out


def human_to_worm(symbols):
    symbols = set(symbols)
    cache = {}
    try:
        cdf = pd.read_csv(ORTH_CACHE, sep="\t")
        for h, sub in cdf.groupby("human_symbol"):
            cache[h] = [w for w in sub["worm_symbol"].dropna().astype(str) if w]
    except FileNotFoundError:
        pass
    missing = [s for s in symbols if s not in cache]
    for i, s in enumerate(missing):
        try:
            cache[s] = _alliance_worm_orthologs(s)
        except Exception:
            cache[s] = []
        time.sleep(0.1)
        if i % 25 == 0:
            print(f"  orthologs {i}/{len(missing)}")
    # rewrite cache
    rows = []
    for h, ws in cache.items():
        if ws:
            rows += [{"human_symbol": h, "worm_symbol": w} for w in ws]
        else:
            rows.append({"human_symbol": h, "worm_symbol": ""})
    pd.DataFrame(rows).sort_values("human_symbol").to_csv(ORTH_CACHE, sep="\t", index=False)
    return {s: cache.get(s, []) for s in symbols}
```

- [ ] **Step 2:** Populate the cache for POS ∪ TOL (long-running first run):
```bash
python -c "import regenerate_grant_figures as g; m=g.human_to_worm(g.load_pos()|g.load_tol()); print('mapped', sum(1 for v in m.values() if v), 'of', len(m))"
```
Expected: prints progress, writes `human_worm_orthologs.tsv`, reports how many human genes have ≥1 worm ortholog.

- [ ] **Step 3:** Verify SON→D1037.1 sanity (known mapping):
```bash
python -c "import pandas as pd; d=pd.read_csv('human_worm_orthologs.tsv',sep='\t'); print(d[d.human_symbol=='SON']); print('rows',len(d))"
```
Expected: SON maps to `D1037.1` (among the worm symbols).

- [ ] **Step 4:** Commit (including the cache):
```bash
git add regenerate_grant_figures.py human_worm_orthologs.tsv && git commit -m "grant figs: human->worm ortholog fetch + cache"
```

---

### Task 3: Truncation-index data + group→TI resolvers

**Files:**
- Modify: `regenerate_grant_figures.py` (add filters, table loaders, resolvers)

**Interfaces:**
- Consumes: `human_to_worm`, `build_worm_groups`, `GRANT_WORM`.
- Produces:
  - `valid_incl0(table) -> DataFrame` (the keep-TI==0 filter).
  - `worm_ti(worm_names: list[str]) -> np.ndarray` — TI values for the transcripts of those worm gene names.
  - `human_ti(symbols: set[str]) -> np.ndarray` — TI values for those human symbols (via ENSG map).
  - `human_all_ti() -> np.ndarray`, `worm_all_ti() -> np.ndarray`.

- [ ] **Step 1:** Add filter + loaders + resolvers:
```python
def valid_incl0(table):
    m = (table["fit_success"] & table["truncationindex"].notna()
         & (table["truncationindex"] >= 0) & (table["truncationindex"] <= 1)
         & (table["n_obs"] >= MIN_OBS))
    return table[m]

_worm_tab = None
_cereb = None  # (valid_table_indexed_by_ENSG, symbol->[ENSG] map)

def _worm_table():
    global _worm_tab
    if _worm_tab is None:
        _worm_tab = pd.read_csv(WORM_TABLE).set_index("gene")
    return _worm_tab

def _cereb_table():
    global _cereb
    if _cereb is None:
        t = pd.read_csv(CEREB_TABLE)
        t["ensg"] = t["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        v = valid_incl0(t)
        src = pd.read_csv(CEREB_CSV, usecols=["Name", "Description"])
        src["ensg"] = src["Name"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        sym2ensg = {}
        for _, r in src.iterrows():
            sym2ensg.setdefault(str(r["Description"]), []).append(r["ensg"])
        _cereb = (v.set_index("ensg"), sym2ensg)
    return _cereb

def worm_all_ti():
    return valid_incl0(_worm_table())["truncationindex"].to_numpy()

def human_all_ti():
    v, _ = _cereb_table()
    return v["truncationindex"].to_numpy()

def worm_ti(worm_names):
    tab = _worm_table()
    sens_like, _ = build_worm_groups_names(worm_names)  # see Step 2
    v = valid_incl0(tab)
    ids = [i for i in sens_like if i in v.index]
    return v.loc[ids, "truncationindex"].to_numpy() if ids else np.array([])

def human_ti(symbols):
    v, sym2ensg = _cereb_table()
    ensgs = [e for s in symbols for e in sym2ensg.get(s, []) if e in v.index]
    ensgs = list(dict.fromkeys(ensgs))
    return v.loc[ensgs, "truncationindex"].to_numpy() if ensgs else np.array([])
```

- [ ] **Step 2:** `build_worm_groups` only resolves the module's `SENS_GENES`/`TOL_GENES`. Add a generalized name→transcript resolver in this file (mirrors the xlsx logic) so any worm-name list works:
```python
def build_worm_groups_names(names):
    """Resolve arbitrary worm gene names -> fourparam transcript IDs present in
    the worm table. Returns (ids, unmapped_names)."""
    from regenerate_acfrog_figures import WORM_MAP_XLSX
    tab = _worm_table()
    xl = pd.read_excel(WORM_MAP_XLSX).dropna(subset=["transcript"])
    xl["tid"] = "w" + xl["wwww"].astype(int).astype(str) + "_" + xl["transcript"].astype(str)
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
        got = by_name.get(n) or by_seq.get(n) or by_seq.get(re.sub(r"\.\d+$", "", str(n)))
        if got:
            ids += got
        else:
            unmapped.append(n)
    return list(dict.fromkeys(ids)), unmapped
```

- [ ] **Step 3:** Verify resolvers return sane sizes:
```bash
python -c "
import regenerate_grant_figures as g
print('worm ALL', len(g.worm_all_ti()), 'human ALL', len(g.human_all_ti()))
print('GRANT worm TI n', len(g.worm_ti(g.GRANT_WORM)))
print('POS human TI n', len(g.human_ti(g.load_pos())))
"
```
Expected: worm ALL in the thousands; human ALL ~15k; GRANT worm TI n > 10; POS human TI n > 20.

- [ ] **Step 4:** Commit:
```bash
git add regenerate_grant_figures.py && git commit -m "grant figs: TI data loaders + group resolvers (keep TI=0)"
```

---

### Task 4: `groups_for(dataset)` assembler

**Files:**
- Modify: `regenerate_grant_figures.py`

**Interfaces:**
- Produces: `groups_for(dataset: str) -> dict[str, np.ndarray]` with keys
  `"POS","GRANT","TOL","ALL"`; `dataset in {"worm","human"}`. Prints a coverage line per group.

- [ ] **Step 1:** Add assembler mapping each set into the dataset's TI space:
```python
def groups_for(dataset):
    pos, tol = load_pos(), load_tol()
    if dataset == "human":
        g = {"POS": human_ti(pos), "GRANT": human_ti(GRANT_HUMAN),
             "TOL": human_ti(tol), "ALL": human_all_ti()}
    elif dataset == "worm":
        h2w = human_to_worm(pos | tol)
        pos_w = [w for s in pos for w in h2w.get(s, [])]
        tol_w = [w for s in tol for w in h2w.get(s, [])]
        g = {"POS": worm_ti(pos_w), "GRANT": worm_ti(GRANT_WORM),
             "TOL": worm_ti(tol_w), "ALL": worm_all_ti()}
    else:
        raise ValueError(dataset)
    for k, v in g.items():
        print(f"  [{dataset}] {k}: n={len(v)}  median={np.median(v):.3f}" if len(v)
              else f"  [{dataset}] {k}: n=0")
    return g
```

- [ ] **Step 2:** Verify both datasets assemble:
```bash
python -c "import regenerate_grant_figures as g; g.groups_for('human'); g.groups_for('worm')"
```
Expected: eight coverage lines, all groups n>0 (POS/GRANT/TOL possibly small in worm).

- [ ] **Step 3:** Commit:
```bash
git add regenerate_grant_figures.py && git commit -m "grant figs: groups_for assembler + coverage report"
```

---

### Task 5: Figure 3B (box/strip + Mann-Whitney), worm + human

**Files:**
- Modify: `regenerate_grant_figures.py`

**Interfaces:**
- Produces: `figure3b(dataset) -> str` (writes `grant_figure3b_<dataset>.png`, returns path).

- [ ] **Step 1:** Implement:
```python
COLORS = {"POS": "#d1495b", "GRANT": "#e3a008", "TOL": "#2e86ab", "ALL": "0.6"}

def _mwu(a, b):
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    return mannwhitneyu(a, b, alternative="two-sided").pvalue

def figure3b(dataset):
    g = groups_for(dataset)
    order = ["POS", "GRANT", "TOL", "ALL"]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    data = [g[k] for k in order]
    bp = ax.boxplot(data, showfliers=False, widths=0.6, patch_artist=True)
    for patch, k in zip(bp["boxes"], order):
        patch.set_facecolor(COLORS[k]); patch.set_alpha(0.35)
    for i, k in enumerate(order, 1):
        y = g[k]
        x = np.random.default_rng(0).normal(i, 0.06, len(y))
        ax.scatter(x, y, s=8, color=COLORS[k], alpha=0.5, linewidths=0)
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels([f"{k}\n(n={len(g[k])})" for k in order])
    ax.set_ylabel("truncation index")
    p_pt = _mwu(g["POS"], g["TOL"]); p_gt = _mwu(g["GRANT"], g["TOL"]); p_pg = _mwu(g["POS"], g["GRANT"])
    ax.set_title(f"Fig 3B ({dataset}): truncation index by OE group\n"
                 f"MWU  POS vs TOL p={p_pt:.3g} | GRANT vs TOL p={p_gt:.3g} | POS vs GRANT p={p_pg:.3g}",
                 fontsize=9)
    fig.tight_layout()
    out = f"grant_figure3b_{dataset}.png"
    fig.savefig(out, dpi=150); plt.close(fig)
    print("wrote", out)
    return out
```

- [ ] **Step 2:** Verify both render:
```bash
python -c "import regenerate_grant_figures as g; g.figure3b('human'); g.figure3b('worm')"
ls -la grant_figure3b_human.png grant_figure3b_worm.png
```
Expected: both PNGs exist, > 20 KB.

- [ ] **Step 3:** Commit:
```bash
git add regenerate_grant_figures.py grant_figure3b_human.png grant_figure3b_worm.png && git commit -m "grant figs: Figure 3B (worm+human)"
```

---

### Task 6: Figure 3D (normalized histograms), worm + human

**Files:**
- Modify: `regenerate_grant_figures.py`

**Interfaces:**
- Produces: `figure3d(dataset) -> str` (writes `grant_figure3d_<dataset>.png`).

- [ ] **Step 1:** Implement per-bin-fraction step histograms of POS/GRANT/TOL vs ALL:
```python
def _prob_hist(ax, vals, color, label, bins):
    if len(vals) == 0:
        return
    counts, edges = np.histogram(vals, bins=bins, range=(0, 1))
    frac = counts / counts.sum()
    centers = (edges[:-1] + edges[1:]) / 2
    ax.step(centers, frac, where="mid", color=color, label=f"{label} (n={len(vals)})", lw=2)

def figure3d(dataset):
    g = groups_for(dataset)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    _prob_hist(ax, g["ALL"], COLORS["ALL"], "All genes", 20)
    _prob_hist(ax, g["TOL"], COLORS["TOL"], "TOL (dup-tolerant)", 20)
    _prob_hist(ax, g["GRANT"], COLORS["GRANT"], "GRANT (Fig 2A)", 20)
    _prob_hist(ax, g["POS"], COLORS["POS"], "POS (positives)", 20)
    ax.set_xlabel("truncation index"); ax.set_ylabel("fraction of genes (p)")
    ax.set_title(f"Fig 3D ({dataset}): normalized truncation-index histograms", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = f"grant_figure3d_{dataset}.png"
    fig.savefig(out, dpi=150); plt.close(fig)
    print("wrote", out)
    return out
```

- [ ] **Step 2:** Verify:
```bash
python -c "import regenerate_grant_figures as g; g.figure3d('human'); g.figure3d('worm')"
ls -la grant_figure3d_human.png grant_figure3d_worm.png
```
Expected: both PNGs exist, > 20 KB.

- [ ] **Step 3:** Commit:
```bash
git add regenerate_grant_figures.py grant_figure3d_human.png grant_figure3d_worm.png && git commit -m "grant figs: Figure 3D (worm+human)"
```

---

### Task 7: Figure 10A (cumulative CDF + KS), human

**Files:**
- Modify: `regenerate_grant_figures.py`

**Interfaces:**
- Produces: `figure10a() -> str` (writes `grant_figure10a_human.png`).

- [ ] **Step 1:** Implement CDF curves + KS vs ALL + uncapped fraction:
```python
def _cdf(ax, vals, color, label):
    if len(vals) == 0:
        return
    x = np.sort(vals)
    y = np.arange(1, len(x) + 1) / len(x)
    frac0 = float(np.mean(vals == 0))
    ax.plot(x, y, color=color, lw=2, label=f"{label} (n={len(vals)}, frac0={frac0:.2f})")

def figure10a():
    g = groups_for("human")
    fig, ax = plt.subplots(figsize=(6.5, 5))
    _cdf(ax, g["ALL"], "k", "All cerebellar genes")
    _cdf(ax, g["TOL"], COLORS["TOL"], "TOL (dup-tolerant)")
    _cdf(ax, g["GRANT"], COLORS["GRANT"], "GRANT (Fig 2A)")
    _cdf(ax, g["POS"], COLORS["POS"], "POS (positives)")
    def ks(a):
        if len(a) < 3:
            return float("nan")
        return ks_2samp(a, g["ALL"]).pvalue
    ax.set_xlabel("truncation index"); ax.set_ylabel("cumulative fraction")
    ax.set_title("Fig 10A (human cerebellum): cumulative truncation index\n"
                 f"KS vs ALL  POS p={ks(g['POS']):.3g} | GRANT p={ks(g['GRANT']):.3g} | TOL p={ks(g['TOL']):.3g}",
                 fontsize=9)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    out = "grant_figure10a_human.png"
    fig.savefig(out, dpi=150); plt.close(fig)
    print("wrote", out)
    return out
```

- [ ] **Step 2:** Verify:
```bash
python -c "import regenerate_grant_figures as g; g.figure10a()"
ls -la grant_figure10a_human.png
```
Expected: PNG exists, > 20 KB; TOL frac0 should be the highest (most uncapped) if the grant pattern holds.

- [ ] **Step 3:** Commit:
```bash
git add regenerate_grant_figures.py grant_figure10a_human.png && git commit -m "grant figs: Figure 10A (human CDF)"
```

---

### Task 8: `__main__` + docs + final run

**Files:**
- Modify: `regenerate_grant_figures.py` (add `__main__`), `CLAUDE.md` (document the new script/figures)

- [ ] **Step 1:** Add runner:
```python
if __name__ == "__main__":
    for ds in ("human", "worm"):
        figure3b(ds); figure3d(ds)
    figure10a()
    print("done: 5 figures")
```

- [ ] **Step 2:** Full run from clean:
```bash
python regenerate_grant_figures.py
ls -la grant_figure3b_worm.png grant_figure3b_human.png grant_figure3d_worm.png grant_figure3d_human.png grant_figure10a_human.png
```
Expected: five PNGs, coverage + stats printed.

- [ ] **Step 3:** Document in `CLAUDE.md` (new section after the grant-proposal section): the script name, the five outputs, the four gene sets, the keep-TI=0 filter divergence, the `human_worm_orthologs.tsv` cache, and the empirical result (fill in actual p-values from the run).

- [ ] **Step 4:** Commit + push:
```bash
git add regenerate_grant_figures.py CLAUDE.md && git commit -m "grant figs: __main__ runner + CLAUDE.md docs" && git push origin main
```

## Self-Review notes

- Spec coverage: POS/GRANT/TOL/ALL (Task 1,4) ✓; ortholog bridge+cache (Task 2) ✓; keep-TI=0 filter (Task 3) ✓; ENSG→symbol map (Task 3) ✓; 3B/3D worm+human + 10A human (Tasks 5-7) ✓; new script, originals untouched ✓.
- The Alliance search-by-symbol endpoint shape is assumed; Task 2 Step 2/3 verify it and the SON→D1037.1 sanity check catches a wrong endpoint early. If the `/api/search` shape differs, fall back to the proven `xrefs`+`/orthologs` path used earlier in the session (resolve human symbol→HGNC/Ensembl id via `rest.ensembl.org/xrefs/symbol/homo_sapiens/<sym>`, then Alliance `/orthologs`, filter taxon C. elegans).
