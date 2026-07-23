#!/usr/bin/env python
"""
generate_positive_controls_fourparam.py

Small, focused fourparam tables for the **key positive-control genes**
(APP, SNCA, PCSK9, SOX9) -- the same 4-parameter Gaussian fit + full column set
as the genome-wide generators, but restricted to these 4 genes so they can be
inspected on their own.

Six tables are written (2 filters x 3 datasets):

    2 filters   : unfiltered (keep every finite value) and excluded <= -1
    3 datasets  : worm, GTEx v10 cerebellum, GTEx v8 cerebellum

    positive_controls_worm_fourparam_table.csv
    positive_controls_worm_fourparam_table_excluded_at_or_below_-1.csv
    positive_controls_cerebellumlog2_fourparam_table.csv
    positive_controls_cerebellumlog2_fourparam_table_excluded_at_or_below_-1.csv
    positive_controls_cerebellumlog2_v8_fourparam_table.csv
    positive_controls_cerebellumlog2_v8_fourparam_table_excluded_at_or_below_-1.csv

Plus two **raw per-sample expression** tables (cerebellum only -- worm lacks an
ortholog for 3 of the 4 controls): the exact rows for the control ENSGs pulled
verbatim from each source matrix, genes-as-rows, already log2-transformed, all
finite samples kept (no exclusion), with a leading ``symbol`` column:

    positive_controls_cerebellumlog2_expression.csv       (4 genes x 266 samples)
    positive_controls_cerebellumlog2_v8_expression.csv    (4 genes x 241 samples)

Each fourparam row is prefixed with a ``symbol`` column (the human positive-control gene)
followed by the exact columns ``BhuvanFitter.fit("fourparam")`` returns (imported
``COLUMNS`` -- so ``gene`` holds the dataset-specific id: an ENSG for cerebellum,
a worm transcript id for worm).

Gene-identifier resolution
--------------------------
* Cerebellum (both versions): human symbol -> versioned ENSG, read straight from
  each cerebellum CSV's own ``Name``/``Description`` columns.
* Worm: human symbol -> *C. elegans* ortholog -> transcript ids. Only **APP ->
  apl-1** has an ortholog in this dataset (SNCA/SOX9 return no worm ortholog from
  the Alliance/DIOPT cache; PCSK9 has no clean 1:1 worm ortholog). The three
  unresolved symbols are emitted as clearly-marked ``fit_success=False`` rows so
  the worm table still lists all four controls.

The fit logic mirrors both genome-wide generators exactly: values are the finite
per-sample expression values (for the excluded tables, values ``<= -1`` are
dropped first); genes with ``< MIN_OBS`` observations or a non-converging fit
get a ``_failed_row``. No refit shortcuts -- every row is a genuine curve_fit.

Run: python generate_positive_controls_fourparam.py [--no-push]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import subprocess

from bhuvanfitter import BhuvanFitter
from generate_fourparam_stats import COLUMNS, MIN_OBS, _failed_row

HERE = Path(__file__).resolve().parent
TABLES = HERE / "outputs/positive_controls"
MAX_NFEV = 2000  # same cap used for the genome-wide cerebellum runs

# Human positive-control symbols, in the order the user listed them.
CONTROLS = ["APP", "SNCA", "PCSK9", "SOX9"]

# ---- per-dataset gene-identifier maps ---------------------------------------
# Worm: only APP has a C. elegans ortholog (apl-1) present in the dataset.
WORM_IDS = {
    "APP": ["w226_C42D8.8a.1", "w227_C42D8.8a.2", "w228_C42D8.8b.1"],  # apl-1
    "SNCA": [],   # synuclein is vertebrate-specific: no C. elegans ortholog
    "PCSK9": [],  # no clean 1:1 worm ortholog (not in Alliance/DIOPT cache)
    "SOX9": [],   # SoxE group: no C. elegans ortholog (worm sox-2/3/4 are other subfamilies)
}
CEREB_V10_IDS = {
    "APP": "ENSG00000142192.22",
    "SNCA": "ENSG00000145335.17",
    "PCSK9": "ENSG00000169174.13",
    "SOX9": "ENSG00000125398.8",
}
CEREB_V8_IDS = {
    "APP": "ENSG00000142192.20",
    "SNCA": "ENSG00000145335.15",
    "PCSK9": "ENSG00000169174.10",
    "SOX9": "ENSG00000125398.5",
}

# (dataset stem, source CSV, leading non-sample columns, symbol->id map)
DATASETS = [
    ("worm", HERE / "data/worm.csv", 1, WORM_IDS),
    ("cerebellumlog2", HERE / "data/cerebellumlog2.csv", 2, CEREB_V10_IDS),
    ("cerebellumlog2_v8", HERE / "data/cerebellumlog2_v8.csv", 2, CEREB_V8_IDS),
]
# (threshold, filename suffix)
THRESHOLDS = [(None, ""), (-1.0, "_excluded_at_or_below_-1")]

OUT_COLUMNS = ["symbol"] + COLUMNS


def fetch_values(csv_path: Path, id_set: set[str], n_lead: int) -> dict[str, np.ndarray]:
    """Single streaming pass over a big genes-as-rows CSV: return the finite-safe
    numeric sample values for every gene id in ``id_set`` (first field == id).
    ``n_lead`` leading non-sample columns are skipped."""
    found: dict[str, np.ndarray] = {}
    remaining = set(id_set)
    if not remaining:
        return found
    with open(csv_path) as f:
        f.readline()  # header
        for line in f:
            comma = line.find(",")
            gid = line[:comma]
            if gid in remaining:
                parts = line.rstrip("\n").split(",")
                found[gid] = np.array([float(x) for x in parts[n_lead:]], dtype=float)
                remaining.discard(gid)
                if not remaining:
                    break
    if remaining:
        raise KeyError(f"{sorted(remaining)} not found in {csv_path.name}")
    return found


def raw_header_and_rows(csv_path: Path, id_set: set[str], n_lead: int):
    """Stream the source matrix once: return (lead_cols, sample_cols, {gid: [str
    fields incl. lead cols]}) so the exact original per-sample values can be
    re-emitted verbatim (no float round-trip)."""
    rows: dict[str, list[str]] = {}
    remaining = set(id_set)
    with open(csv_path) as f:
        header = f.readline().rstrip("\n").split(",")
        lead_cols, sample_cols = header[:n_lead], header[n_lead:]
        if remaining:
            for line in f:
                comma = line.find(",")
                gid = line[:comma]
                if gid in remaining:
                    rows[gid] = line.rstrip("\n").split(",")
                    remaining.discard(gid)
                    if not remaining:
                        break
    if remaining:
        raise KeyError(f"{sorted(remaining)} not found in {csv_path.name}")
    return lead_cols, sample_cols, rows


def write_raw_expression(name: str, csv_path: Path, n_lead: int, id_map: dict) -> Path:
    """Write the raw (already log-transformed) per-sample expression rows for the
    control genes, genes-as-rows, mirroring the source matrix with a leading
    ``symbol`` column. All finite values kept -- this is the unfiltered raw data."""
    ids = {gid for gid in id_map.values() if isinstance(gid, str) and gid}
    lead_cols, sample_cols, src = raw_header_and_rows(csv_path, ids, n_lead)
    records = []
    for symbol in CONTROLS:
        gid = id_map[symbol]
        if not (isinstance(gid, str) and gid):
            continue
        fields = src[gid]
        rec = {"symbol": symbol}
        rec.update(zip(lead_cols, fields[:n_lead]))
        rec.update(zip(sample_cols, fields[n_lead:]))
        records.append(rec)
    out_cols = ["symbol"] + lead_cols + sample_cols
    table = pd.DataFrame.from_records(records, columns=out_cols)
    out = TABLES / f"positive_controls_{name}_expression.csv"
    table.to_csv(out, index=False)
    print(f"  wrote {out.name}  ({len(table)} genes x {len(sample_cols)} samples)")
    return out


def fit_values(gid: str, values: np.ndarray, threshold) -> dict:
    """Fit one gene's value array (mirrors both generators' _fit_one)."""
    data = values[np.isfinite(values)]
    if threshold is not None:
        data = data[data > threshold]  # drop <= threshold, same as the excluded generator
    n_obs = int(data.size)
    if n_obs < MIN_OBS:
        return _failed_row(gid, n_obs)
    try:
        bf = BhuvanFitter(data, gene_name=gid)
        return bf.fit("fourparam", max_nfev=MAX_NFEV)
    except RuntimeError:
        return _failed_row(gid, n_obs)


def no_ortholog_row(symbol: str) -> dict:
    """A clearly-marked placeholder row for a control with no id in this dataset."""
    row = _failed_row("(no worm ortholog)", 0)
    row["symbol"] = symbol
    return row


def build_dataset_tables(name: str, csv_path: Path, n_lead: int, id_map: dict,
                         no_push: bool) -> list[Path]:
    """Write the unfiltered + excluded tables for one dataset."""
    # Collect every id we need for this dataset in one streaming pass.
    all_ids: set[str] = set()
    for gid in id_map.values():
        if isinstance(gid, list):
            all_ids.update(gid)
        elif gid:
            all_ids.add(gid)
    print(f"\n=== {name}: reading {len(all_ids)} gene rows from {csv_path.name} ===")
    values = fetch_values(csv_path, all_ids, n_lead)

    written: list[Path] = []
    for threshold, suffix in THRESHOLDS:
        rows: list[dict] = []
        for symbol in CONTROLS:
            gid = id_map[symbol]
            ids = gid if isinstance(gid, list) else ([gid] if gid else [])
            if not ids:
                rows.append(no_ortholog_row(symbol))
                continue
            for one in ids:
                row = fit_values(one, values[one], threshold)
                row = {"symbol": symbol, **row}
                rows.append(row)

        table = pd.DataFrame.from_records(rows, columns=OUT_COLUMNS)
        out = TABLES / f"positive_controls_{name}_fourparam_table{suffix}.csv"
        table.to_csv(out, index=False)
        label = "unfiltered" if threshold is None else f"excluded <= {threshold:g}"
        n_ok = int(table["fit_success"].sum())
        print(f"  wrote {out.name}  ({label}: {n_ok}/{len(table)} rows fit)")
        written.append(out)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-push", action="store_true",
                        help="Write the tables but do not commit/push to GitHub.")
    args = parser.parse_args()

    TABLES.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, csv_path, n_lead, id_map in DATASETS:
        written.extend(build_dataset_tables(name, csv_path, n_lead, id_map,
                                            args.no_push))

    # Raw per-sample (already log-transformed) expression for the control genes,
    # cerebellum datasets only (worm has no ortholog for 3 of the 4 controls).
    print("\n=== raw expression tables (cerebellum) ===")
    for name, csv_path, n_lead, id_map in DATASETS:
        if not name.startswith("cerebellum"):
            continue
        written.append(write_raw_expression(name, csv_path, n_lead, id_map))

    # Compact TI summary across every fourparam table written (raw expression
    # tables have no truncationindex column, so they're skipped).
    print("\nTruncation index by table:")
    for out in written:
        if "fourparam" not in out.name:
            continue
        t = pd.read_csv(out)
        for _, r in t.iterrows():
            ti = r["truncationindex"]
            ti_s = f"{ti:.4f}" if pd.notna(ti) else "  n/a "
            print(f"  {out.name:66s}  {str(r['symbol']):6s} "
                  f"{str(r['gene']):22s} TI={ti_s}")

    if args.no_push:
        print("\n--no-push set; skipping git push.")
        return
    push_tables(written)


def push_tables(paths: list[Path]) -> None:
    """Stage the given tables (by their real repo-relative paths -- they live in a
    subfolder, so a basename-only ``git add`` would miss them), commit, and push.
    Only these files are staged, so unrelated working-tree changes are untouched."""
    def run(*args: str) -> str:
        r = subprocess.run(["git", "-C", str(HERE), *args],
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"`git {' '.join(args)}` failed:\n{r.stderr.strip()}")
        return r.stdout

    rel = [str(p.relative_to(HERE).as_posix()) for p in paths]
    run("add", *rel)
    if not run("status", "--porcelain", *rel).strip():
        print("No table changes to commit or push.")
        return
    run("commit", "-m", "Update positive-control fourparam tables")
    run("push", "origin", "HEAD")
    print(f"Pushed {len(rel)} positive-control tables to origin.")


if __name__ == "__main__":
    main()
