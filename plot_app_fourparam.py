#!/usr/bin/env python
"""
plot_app_fourparam.py

Histogram + fitted 4-parameter Gaussian curve for the APP gene across every
worm and cerebellum fourparam table, using the **already-generated** y0/A/x0/w
from each table (no refit). The truncation index (from the same table row) is
printed under each panel.

APP is human; its C. elegans ortholog is apl-1 (transcripts C42D8.8a.1 /
C42D8.8a.2 / C42D8.8b.1). The fourparam curve was fit to the 40-bin histogram
*counts*, so ``f(x) = y0 + A*exp(-((x-x0)/w)^2)`` is drawn directly on the count
axis (matching BhuvanFitter.hist), over each table's exact input array (finite
values, with the table's ``<= threshold`` exclusion applied).

Run: python plot_app_fourparam.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TABLES = HERE / "outputs/tables"
FIGS = HERE / "outputs/figures"
BINS = 40


def fetch_row_values(csv_path: Path, gene_id: str, n_lead: int) -> np.ndarray:
    """Stream a big genes-as-rows CSV and return the numeric sample values for
    one gene (first field == gene_id). ``n_lead`` non-sample leading columns are
    skipped (worm: 1 = strain; GTEx: 2 = Name,Description)."""
    with open(csv_path) as f:
        f.readline()  # header
        for line in f:
            comma = line.find(",")
            if line[:comma] == gene_id:
                parts = line.rstrip("\n").split(",")
                return np.array([float(x) for x in parts[n_lead:]], dtype=float)
    raise KeyError(f"{gene_id} not found in {csv_path.name}")


def table_params(table_csv: Path, gene_id: str) -> dict:
    """Pull the stored fit row (y0/A/x0/w/truncationindex/n_obs) for one gene."""
    t = pd.read_csv(table_csv)
    r = t[t["gene"].astype(str) == gene_id]
    if r.empty:
        raise KeyError(f"{gene_id} not in {table_csv.name}")
    return r.iloc[0].to_dict()


def curve(x, p) -> np.ndarray:
    return p["y0"] + p["A"] * np.exp(-((x - p["x0"]) / p["w"]) ** 2)


def draw_panel(ax, values: np.ndarray, threshold, p: dict, title: str) -> None:
    data = values[np.isfinite(values)]
    if threshold is not None:
        data = data[data > threshold]  # same exclusion the table used
    counts, edges = np.histogram(data, bins=BINS)
    centers = (edges[:-1] + edges[1:]) / 2
    ax.bar(centers, counts, width=np.diff(edges).mean(),
           color="steelblue", alpha=0.6, label="Observed", zorder=2)
    xs = np.linspace(edges[0], edges[-1], 600)
    ax.plot(xs, curve(xs, p), color="crimson", linewidth=2,
            label="4-param Gaussian", zorder=3)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(f"expression   (truncation index = {p['truncationindex']:.4f})",
                  fontsize=9)
    ax.set_ylabel("count", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")


# ---- WORM: apl-1 (APP ortholog), 3 transcripts x 3 tables --------------------
WORM_CSV = HERE / "data/worm.csv"
WORM_TX = ["w226_C42D8.8a.1", "w227_C42D8.8a.2", "w228_C42D8.8b.1"]
WORM_TABLES = [
    ("worm_fourparam_table.csv", "unfiltered", None),
    ("worm_fourparam_table_excluded_at_or_below_-0.75.csv", "excl ≤ −0.75", -0.75),
    ("worm_fourparam_table_excluded_at_or_below_-1.csv", "excl ≤ −1", -1.0),
]


def plot_worm() -> Path:
    vals = {tx: fetch_row_values(WORM_CSV, tx, n_lead=1) for tx in WORM_TX}
    fig, axes = plt.subplots(len(WORM_TABLES), len(WORM_TX),
                             figsize=(15, 11), squeeze=False)
    for i, (fname, tlabel, thr) in enumerate(WORM_TABLES):
        for j, tx in enumerate(WORM_TX):
            p = table_params(TABLES / fname, tx)
            title = f"{tx}  (apl-1)\n{tlabel}"
            draw_panel(axes[i][j], vals[tx], thr, p, title)
    fig.suptitle(
        "APP worm ortholog apl-1 — histogram + generated 4-param Gaussian, all 3 worm tables\n"
        "(rows are identical: apl-1 has no expression values ≤ −0.75, so exclusion removes nothing; n_obs=207)",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIGS / "app_apl1_worm_fourparam.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# ---- CEREBELLUM: APP, 3 tables ----------------------------------------------
CEREB_PANELS = [
    ("cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv", "GTEx v10  —  excl ≤ −1",
     -1.0, HERE / "data/cerebellumlog2.csv", "ENSG00000142192.22"),
    ("cerebellumlog2_v8_fourparam_table.csv", "GTEx v8  —  unfiltered",
     None, HERE / "data/cerebellumlog2_v8.csv", "ENSG00000142192.20"),
    ("cerebellumlog2_v8_fourparam_table_excluded_at_or_below_-1.csv", "GTEx v8  —  excl ≤ −1",
     -1.0, HERE / "data/cerebellumlog2_v8.csv", "ENSG00000142192.20"),
]


def plot_cerebellum() -> Path:
    fig, axes = plt.subplots(1, len(CEREB_PANELS), figsize=(16, 5), squeeze=False)
    for j, (fname, tlabel, thr, src, gid) in enumerate(CEREB_PANELS):
        vals = fetch_row_values(src, gid, n_lead=2)
        p = table_params(TABLES / fname, gid)
        title = f"APP  ({gid})\n{tlabel}"
        draw_panel(axes[0][j], vals, thr, p, title)
    fig.suptitle(
        "APP (ENSG00000142192) — histogram + generated 4-param Gaussian, all cerebellum tables\n"
        "(APP is highly expressed: min log2 value well above −1, so the v8 unfiltered and v8 excl≤−1 fits are identical)",
        fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out = FIGS / "app_cerebellum_fourparam.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    w = plot_worm()
    c = plot_cerebellum()
    print(f"Wrote {w}")
    print(f"Wrote {c}")

    # Print the TI table for all panels.
    print("\nTruncation index (from the generated tables):")
    for fname, tlabel, _ in WORM_TABLES:
        for tx in WORM_TX:
            p = table_params(TABLES / fname, tx)
            print(f"  worm  {tlabel:14s}  {tx:18s} apl-1   TI = {p['truncationindex']:.4f}")
    for fname, tlabel, _, _, gid in CEREB_PANELS:
        p = table_params(TABLES / fname, gid)
        print(f"  cereb {tlabel:22s}  {gid:20s} APP   TI = {p['truncationindex']:.4f}")


if __name__ == "__main__":
    main()
