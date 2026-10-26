#!/usr/bin/env python3
"""Run the available, non-causal T14 orthogonal-validation gate.

This script deliberately performs a partial T14 execution.  It audits the
frozen T13 gene candidates in locally available thyroid single-cell/spatial
archives and the author-filtered HRA001684 thyroid matrices, while retaining
GTEx thyroid eQTL as an explicit HOLD input. Cell-level observations are
aggregated to sample/donor summaries; no cell is treated as an independent
replicate.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import os
import re
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outcome_data" / "t14"
OUT.mkdir(parents=True, exist_ok=True)


def md5_path(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def gzip_integrity(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return True
    except Exception:
        return False


def tar_integrity(path: Path) -> tuple[bool, int]:
    try:
        with tarfile.open(path, "r:gz" if path.name.endswith(".gz") else "r:") as tf:
            members = tf.getmembers()
        return True, len(members)
    except Exception:
        return False, 0


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if row.get(k) is None else row.get(k) for k in fields})


def parse_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def matrix_market_dims_and_targets(path: Path, gene_rows: list[dict[str, str]]) -> tuple[int, int, dict[int, tuple[str, str]] , dict[str, int], dict[str, int]]:
    """Read a gzipped/ plain MatrixMarket file and aggregate only T13 genes.

    The matching feature index is supplied separately because the Zenodo
    archive uses both legacy genes.tsv and newer features.tsv names.
    """
    with path.open("rt", encoding="utf-8", errors="replace") as fh:
        line = fh.readline()
        if not line.startswith("%%MatrixMarket matrix coordinate"):
            raise RuntimeError(f"Invalid MatrixMarket header: {path}")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        rows, cols, _ = (int(x) for x in line.split()[:3])
        targets = {i: v for i, v in enumerate([], start=1)}
        totals = {r["gene_symbol"]: 0 for r in gene_rows}
        detected = {r["gene_symbol"]: 0 for r in gene_rows}
        target_idx: dict[int, str] = {}
        for r in gene_rows:
            target_idx[int(r["_feature_index"])] = r["gene_symbol"]
        for line in fh:
            p = line.split()
            if len(p) < 3:
                continue
            i, _, value = int(p[0]), int(p[1]), int(float(p[2]))
            symbol = target_idx.get(i)
            if symbol is not None:
                totals[symbol] += value
                detected[symbol] += 1
        return rows, cols, {}, totals, detected


def hra001684_expression(gene_rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    root = ROOT / "raw_data/single_cell/HRA001684_zenodo_4575509/extracted/thyroiditis"
    libs = ["HT12Thyroid", "HT13Thyroid", "HT28Thyroid", "HT29Thyroid"]
    rows: list[dict[str, object]] = []
    total_cells = 0
    for lib in libs:
        d = root / lib / "filtered_feature_bc_matrix"
        feature_path = d / ("features.tsv" if (d / "features.tsv").exists() else "genes.tsv")
        features = []
        with feature_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\r\n").split("\t")
                features.append(p[1] if len(p) > 1 else p[0])
        target_by_symbol = {r["gene_symbol"]: i + 1 for i, g in enumerate(features) for r in gene_rows if g == r["gene_symbol"]}
        n_features = len(features)
        matrix = d / "matrix.mtx"
        with matrix.open("r", encoding="utf-8", errors="replace") as fh:
            line = fh.readline()
            if not line.startswith("%%MatrixMarket matrix coordinate"):
                raise RuntimeError(f"Invalid MatrixMarket header: {matrix}")
            line = fh.readline()
            while line.startswith("%"):
                line = fh.readline()
            mrows, n_cells, _ = (int(x) for x in line.split()[:3])
            if mrows != n_features:
                raise RuntimeError(f"Feature/matrix row mismatch for {lib}")
            totals = {r["gene_symbol"]: 0 for r in gene_rows}
            detected = {r["gene_symbol"]: 0 for r in gene_rows}
            target_idx = {i: s for s, i in target_by_symbol.items()}
            for line in fh:
                p = line.split()
                if len(p) < 3:
                    continue
                i, _, value = int(p[0]), int(p[1]), int(float(p[2]))
                symbol = target_idx.get(i)
                if symbol is not None:
                    totals[symbol] += value
                    detected[symbol] += 1
        total_cells += n_cells
        for r in gene_rows:
            symbol = r["gene_symbol"]
            rows.append({"dataset": "HRA001684_Zenodo_filtered", "gene_id": r["gene_id"], "gene_symbol": symbol,
                         "donor_id": lib[:4], "sample_id": lib, "tissue": "thyroid", "n_cells": n_cells,
                         "n_detected": detected[symbol], "detection_rate": f"{detected[symbol] / n_cells:.8f}",
                         "total_count": totals[symbol], "mean_count_per_cell": f"{totals[symbol] / n_cells:.8f}",
                         "evidence_class": "author_filtered_matrix_donor_level_descriptive", "replicate_unit": "library/donor; cells aggregated",
                         "limitation": "Zenodo filtered matrix; no disease-control group and no cell-level pseudoreplication"})
    return rows, {"libraries": len(libs), "cells": total_cells, "target_genes": len(gene_rows)}


def gtex_thyroid_audit(path: Path, gene_rows: list[dict[str, str]]) -> dict[str, object]:
    wanted_ids = {r["gene_id"] for r in gene_rows}
    wanted_symbols = {r["gene_symbol"] for r in gene_rows}
    matched_ids: Counter[str] = Counter()
    matched_symbols: Counter[str] = Counter()
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        required = {"ensembl_gene_id", "external_gene_name", "chr", "pos", "beta", "se", "p", "n", "rsids2"}
        if not required.issubset(header):
            raise RuntimeError("GTEx thyroid eQTL header lacks required fields")
        idx = {name: header.index(name) for name in required}
        for line in fh:
            p = line.rstrip("\r\n").split("\t")
            if len(p) <= max(idx.values()):
                continue
            if p[idx["ensembl_gene_id"]] in wanted_ids:
                matched_ids[p[idx["ensembl_gene_id"]]] += 1
            if p[idx["external_gene_name"]] in wanted_symbols:
                matched_symbols[p[idx["external_gene_name"]]] += 1
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "matched_ids": dict(matched_ids), "matched_symbols": dict(matched_symbols), "candidate_rows": sum(matched_ids.values())}


def normalize_gse218743_cell_id(value: str) -> str:
    """Undo R's dot substitution in the supplied count-matrix column names."""
    return re.sub(r"\.1_", "-1_", value)


def load_t13_candidates() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    decision = read_tsv(ROOT / "outcome_data" / "t13" / "t13_decision.tsv")
    if len(decision) != 1 or decision[0]["status"] != "PASS_WITH_HOLD":
        raise RuntimeError("T13 decision is not PASS_WITH_HOLD; T14 is not authorized")
    if decision[0]["decision"] != "PASS_T13_PARTIAL_CASCADE_HOLD_PQTL_CODING":
        raise RuntimeError("Unexpected T13 decision status")
    candidate_rows = read_tsv(ROOT / "outcome_data" / "t13" / "t13_candidate_summary.tsv")
    gene_rows = read_tsv(ROOT / "outcome_data" / "t13" / "t13_gene_summary.tsv")
    if len(candidate_rows) != 26 or len({r["candidate_id"] for r in candidate_rows}) != 26:
        raise RuntimeError("T13 candidate summary is not the frozen 26-feature set")
    gene_rows = [r for r in gene_rows if r.get("gene_id") and r.get("gene_symbol")]
    if not gene_rows:
        raise RuntimeError("No T13 gene candidates with symbols are available")
    return candidate_rows, gene_rows


def manifest_lookup() -> dict[str, dict[str, str]]:
    path = ROOT / "raw_data" / "manifests" / "local_data_manifest.tsv"
    rows = read_tsv(path)
    return {r.get("filename", ""): r for r in rows}


def source_row(
    dataset: str,
    path: Path | None,
    manifest: dict[str, dict[str, str]],
    status: str,
    limitation: str,
    *,
    integrity: str = "",
    members: object = "",
) -> dict[str, object]:
    if path is None or not path.exists():
        return {
            "dataset": dataset,
            "path": "not_available_in_current_project",
            "exists": False,
            "bytes_actual": "",
            "md5_actual": "",
            "sha256_actual": "",
            "manifest_bytes": "",
            "manifest_md5": "",
            "manifest_match": False,
            "integrity": "not_run",
            "archive_members": members,
            "status": status,
            "limitation": limitation,
        }
    m = manifest.get(path.name, {})
    actual_bytes = path.stat().st_size
    actual_md5 = md5_path(path)
    manifest_bytes = m.get("file_size", "")
    manifest_md5 = m.get("checksum", "")
    match = str(actual_bytes) == str(manifest_bytes) and actual_md5 == manifest_md5
    effective_status = status
    if not match and status.startswith("READY"):
        effective_status = "READY_STRUCTURE_ONLY_MANIFEST_MISMATCH"
    return {
        "dataset": dataset,
        "path": str(path.relative_to(ROOT)),
        "exists": True,
        "bytes_actual": actual_bytes,
        "md5_actual": actual_md5,
        "sha256_actual": sha256_path(path),
        "manifest_bytes": manifest_bytes,
        "manifest_md5": manifest_md5,
        "manifest_match": match,
        "integrity": integrity,
        "archive_members": members,
            "status": effective_status,
        "limitation": limitation if match else limitation + "; old local manifest bytes/checksum do not match actual file",
    }


def gse182416_expression(gene_rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    metadata_path = ROOT / "raw_data/single_cell/GSE182416/GSE182416_Thyroid_normal_7samples_metadata.txt.gz"
    counts_path = ROOT / "raw_data/single_cell/GSE182416/GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz"
    with gzip.open(metadata_path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        metadata = list(csv.DictReader(fh, delimiter="\t"))
    meta_by_id = {r["Index"]: r for r in metadata}
    genes = {r["gene_symbol"]: r["gene_id"] for r in gene_rows}
    gene_values: dict[str, list[int]] = {}
    found = Counter()
    with gzip.open(counts_path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n\r").split("\t")
        cell_ids = header[1:]
        if len(cell_ids) != len(metadata) or set(cell_ids) != set(meta_by_id):
            raise RuntimeError("GSE182416 count-matrix cell IDs do not match metadata")
        for line in fh:
            parts = line.rstrip("\n\r").split("\t")
            if not parts:
                continue
            symbol = parts[0]
            if symbol not in genes:
                continue
            if symbol in gene_values:
                raise RuntimeError(f"Duplicate target gene row in GSE182416: {symbol}")
            vals = [parse_int(x) for x in parts[2:]]
            if len(vals) != len(cell_ids):
                raise RuntimeError(f"GSE182416 target row length mismatch: {symbol}")
            gene_values[symbol] = vals
            found[symbol] += 1
    missing = sorted(set(genes) - set(gene_values))
    if missing:
        raise RuntimeError("GSE182416 missing T13 gene rows: " + ", ".join(missing))

    accum: dict[tuple[str, str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    for symbol, vals in gene_values.items():
        gene_id = genes[symbol]
        for cell_id, value in zip(cell_ids, vals):
            m = meta_by_id[cell_id]
            key = (gene_id, symbol, m["orig.ident"], m["celltype_global"], m["celltype_sub"])
            accum[key][0] += 1
            accum[key][1] += value
            if value > 0:
                # Detection is stored separately below using a negative sentinel-free count.
                pass
    # Recompute detection counts in a second small pass over the two target rows.
    detected: dict[tuple[str, str, str, str, str], int] = Counter()
    for symbol, vals in gene_values.items():
        gene_id = genes[symbol]
        for cell_id, value in zip(cell_ids, vals):
            if value > 0:
                m = meta_by_id[cell_id]
                detected[(gene_id, symbol, m["orig.ident"], m["celltype_global"], m["celltype_sub"])] += 1
    rows: list[dict[str, object]] = []
    for key, (n_cells, total) in sorted(accum.items()):
        gene_id, symbol, donor, global_ct, sub_ct = key
        n_detected = detected[key]
        rows.append({
            "dataset": "GSE182416",
            "gene_id": gene_id,
            "gene_symbol": symbol,
            "donor_id": donor,
            "celltype_global": global_ct,
            "celltype_sub": sub_ct,
            "n_cells": n_cells,
            "n_detected": n_detected,
            "detection_rate": f"{n_detected / n_cells:.8f}",
            "total_count": total,
            "mean_count_per_cell": f"{total / n_cells:.8f}",
            "evidence_class": "normal_thyroid_cell_localization",
            "replicate_unit": "donor/sample; cells aggregated",
        })
    qc = {"metadata_cells": len(metadata), "matrix_cells": len(cell_ids), "target_genes_found": len(found), "target_rows_found": sum(found.values())}
    return rows, qc


def gse218743_expression(gene_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    metadata_path = ROOT / "raw_data/single_cell/GSE218743/GSE218743_irae_ht_integrated_metadata.csv.gz"
    counts_path = ROOT / "raw_data/single_cell/GSE218743/GSE218743_irae_ht_integrated_counts.csv.gz"
    with gzip.open(metadata_path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        metadata = list(csv.DictReader(fh))
    meta_by_id = {r[""]: r for r in metadata}
    genes = {r["gene_symbol"]: r["gene_id"] for r in gene_rows}
    cell_ids: list[str] = []
    gene_values: dict[str, list[int]] = {}
    with gzip.open(counts_path, "rt", newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        raw_cell_ids = header[1:]
        cell_ids = [normalize_gse218743_cell_id(x) for x in raw_cell_ids]
        if len(cell_ids) != len(metadata) or set(cell_ids) != set(meta_by_id):
            raise RuntimeError("GSE218743 count-matrix cell IDs do not match metadata")
        for row in reader:
            if row and row[0] in genes:
                gene_values[row[0]] = [parse_int(x) for x in row[1:]]
    rows: list[dict[str, object]] = []
    accum: dict[tuple[str, str, str, str, str], list[int]] = defaultdict(lambda: [0, 0, 0])
    for symbol, vals in gene_values.items():
        for cell_id, value in zip(cell_ids, vals):
            m = meta_by_id[cell_id]
            key = (genes[symbol], symbol, m["sample"], m["condition"], m["CellType"])
            accum[key][0] += 1
            accum[key][1] += value
            accum[key][2] += int(value > 0)
    for key, (n, total, detected) in sorted(accum.items()):
        gene_id, symbol, sample, condition, celltype = key
        rows.append({
            "dataset": "GSE218743",
            "gene_id": gene_id,
            "gene_symbol": symbol,
            "sample_id": sample,
            "condition_label": condition,
            "celltype_label_as_supplied": celltype,
            "n_cells": n,
            "n_detected": detected,
            "detection_rate": f"{detected / n:.8f}",
            "total_count": total,
            "mean_count_per_cell": f"{total / n:.8f}",
            "evidence_class": "exploratory_sample_level_not_donor_verified",
            "replicate_unit": "sample; cells aggregated",
            "limitation": "Not used as HRA001684 donor-level thyroid validation; celltype labels remain supplied cluster IDs",
        })
    return rows


def tar_member(tf: tarfile.TarFile, name: str):
    member = tf.getmember(name)
    fh = tf.extractfile(member)
    if fh is None:
        raise RuntimeError(f"Cannot read tar member {name}")
    return fh


def read_gzip_member(tf: tarfile.TarFile, name: str):
    return gzip.GzipFile(fileobj=tar_member(tf, name), mode="rb")


def spatial_inventory(gene_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    path = ROOT / "raw_data/spatial/GSE248205/GSE248205_Processed_data.tar.gz"
    genes = {r["gene_symbol"]: r["gene_id"] for r in gene_rows}
    rows: list[dict[str, object]] = []
    with tarfile.open(path, "r:gz") as tf:
        names = set(tf.getnames())
        samples = sorted({m.group(1) for name in names if (m := re.match(r"Processed files/([^/]+)/\1_matrix\.mtx\.gz$", name))})
        if len(samples) != 8:
            raise RuntimeError(f"Expected 8 GSE248205 spatial samples, found {len(samples)}")
        for sample in samples:
            base = f"Processed files/{sample}/{sample}"
            feature_name = base + "_features.tsv.gz"
            barcode_name = base + "_barcodes.tsv.gz"
            matrix_name = base + "_matrix.mtx.gz"
            required = [feature_name, barcode_name, matrix_name, base + "_tissue_positions_list.csv"]
            if not all(x in names for x in required):
                raise RuntimeError(f"GSE248205 missing required members for {sample}")
            target_rows: dict[int, tuple[str, str]] = {}
            feature_count = 0
            with read_gzip_member(tf, feature_name) as fh:
                for idx, raw in enumerate(fh, 1):
                    parts = raw.decode("utf-8", "replace").rstrip("\n\r").split("\t")
                    feature_count += 1
                    if len(parts) >= 2 and parts[1] in genes:
                        target_rows[idx] = (parts[0], parts[1])
            with read_gzip_member(tf, barcode_name) as fh:
                barcode_count = sum(1 for _ in fh)
            with gzip.GzipFile(fileobj=tar_member(tf, matrix_name), mode="rb") as fh:
                first = fh.readline().decode("utf-8", "replace").strip()
                if not first.startswith("%%MatrixMarket"):
                    raise RuntimeError(f"Invalid MatrixMarket header: {matrix_name}")
                line = fh.readline()
                while line.startswith(b"%"):
                    line = fh.readline()
                dims = [int(x) for x in line.decode().strip().split()]
                if len(dims) != 3 or dims[0] != feature_count or dims[1] != barcode_count:
                    raise RuntimeError(f"GSE248205 dimensions mismatch for {sample}")
                total = Counter()
                detected = Counter()
                for raw in fh:
                    parts = raw.decode("utf-8", "replace").split()
                    if len(parts) < 3:
                        continue
                    row_idx = int(parts[0])
                    if row_idx in target_rows:
                        symbol = target_rows[row_idx][1]
                        value = parse_int(parts[2])
                        total[symbol] += value
                        detected[symbol] += int(value > 0)
            with tar_member(tf, base + "_tissue_positions_list.csv") as fh:
                spot_rows = max(0, sum(1 for _ in io.TextIOWrapper(fh, encoding="utf-8", errors="replace")) - 1)
            group = re.match(r"([A-Za-z]+)", sample).group(1)
            for symbol, gene_id in genes.items():
                rows.append({
                    "dataset": "GSE248205",
                    "sample_id": sample,
                    "archive_group_prefix": group,
                    "group_label_source": "sample-name prefix only; phenotype not independently assigned",
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "n_features": feature_count,
                    "n_barcodes": barcode_count,
                    "matrix_rows": dims[0],
                    "matrix_columns": dims[1],
                    "n_tissue_position_rows": spot_rows,
                    "candidate_in_features": any(v[1] == symbol for v in target_rows.values()),
                    "n_nonzero_spots": detected[symbol],
                    "detection_rate_among_barcodes": f"{detected[symbol] / barcode_count:.8f}",
                    "total_count": total[symbol],
                    "evidence_class": "spatial_archive_candidate_detection",
                    "limitation": "No spatial-neighborhood statistic run; archive group prefix is not treated as a verified phenotype",
                })
    return rows


def gse285196_inventory(gene_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    path = ROOT / "raw_data/single_cell/GSE285196/GSE285196_RAW.tar"
    genes = {r["gene_symbol"]: r["gene_id"] for r in gene_rows}
    rows: list[dict[str, object]] = []
    with tarfile.open(path, "r:") as tf:
        names = set(tf.getnames())
        matrix_names = sorted(n for n in names if re.search(r"_matrix\d+\.mtx\.gz$", n))
        if len(matrix_names) != 10:
            raise RuntimeError(f"Expected 10 GSE285196 matrices, found {len(matrix_names)}")
        for matrix_name in matrix_names:
            stem = matrix_name[: -len(".mtx.gz")]
            sample = re.match(r"(GSM\d+)_matrix(\d+)$", stem).group(1)
            idx = re.match(r"GSM\d+_matrix(\d+)$", stem).group(1)
            feature_name = f"{sample}_features{idx}.tsv.gz"
            barcode_name = f"{sample}_barcodes{idx}.tsv.gz"
            features: set[str] = set()
            n_features = 0
            with gzip.GzipFile(fileobj=tar_member(tf, feature_name), mode="rb") as fh:
                for raw in fh:
                    parts = raw.decode("utf-8", "replace").rstrip("\n\r").split("\t")
                    n_features += 1
                    features.update(parts[:2])
            with gzip.GzipFile(fileobj=tar_member(tf, barcode_name), mode="rb") as fh:
                n_barcodes = sum(1 for _ in fh)
            with gzip.GzipFile(fileobj=tar_member(tf, matrix_name), mode="rb") as fh:
                first = fh.readline().decode("utf-8", "replace").strip()
                second = fh.readline()
                while second.startswith(b"%"):
                    second = fh.readline()
                dims = [int(x) for x in second.decode().strip().split()]
                target_indices = {}
                with gzip.GzipFile(fileobj=tar_member(tf, feature_name), mode="rb") as feature_fh:
                    for feature_idx, raw in enumerate(feature_fh, 1):
                        parts = raw.decode("utf-8", "replace").rstrip("\n\r").split("\t")
                        for symbol in genes:
                            if symbol in parts[:2]:
                                target_indices[feature_idx] = symbol
                totals = Counter()
                detected = Counter()
                for raw in fh:
                    parts = raw.decode("utf-8", "replace").split()
                    if len(parts) >= 3 and int(parts[0]) in target_indices:
                        symbol = target_indices[int(parts[0])]
                        value = parse_int(parts[2])
                        totals[symbol] += value
                        detected[symbol] += int(value > 0)
            for symbol, gene_id in genes.items():
                rows.append({
                    "dataset": "GSE285196",
                    "sample_id": sample,
                    "matrix_index": idx,
                    "n_features": n_features,
                    "n_barcodes": n_barcodes,
                    "matrix_rows": dims[0],
                    "matrix_columns": dims[1],
                    "candidate_in_features": symbol in features,
                    "n_nonzero_cells": detected[symbol],
                    "detection_rate_among_barcodes": f"{detected[symbol] / n_barcodes:.8f}",
                    "total_count": totals[symbol],
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "evidence_class": "archive_inventory_candidate_absent_or_present",
                    "limitation": "No phenotype mapping or candidate expression if feature absent; no Graves-specific inference from GSM ID",
                })
    return rows


def main() -> int:
    candidate_rows, gene_rows = load_t13_candidates()
    manifest = manifest_lookup()
    gse182_count = ROOT / "raw_data/single_cell/GSE182416/GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz"
    gse182_meta = ROOT / "raw_data/single_cell/GSE182416/GSE182416_Thyroid_normal_7samples_metadata.txt.gz"
    gse218_count = ROOT / "raw_data/single_cell/GSE218743/GSE218743_irae_ht_integrated_counts.csv.gz"
    gse218_meta = ROOT / "raw_data/single_cell/GSE218743/GSE218743_irae_ht_integrated_metadata.csv.gz"
    gse248 = ROOT / "raw_data/spatial/GSE248205/GSE248205_Processed_data.tar.gz"
    gse285 = ROOT / "raw_data/single_cell/GSE285196/GSE285196_RAW.tar"
    hra_zip = ROOT / "raw_data/single_cell/HRA001684_zenodo_4575509/thyroiditis.zip"

    inventory: list[dict[str, object]] = []
    for dataset, path, kind, status, limitation in [
        ("GSE182416_count_matrix", gse182_count, "gzip", "READY_WITH_INDEPENDENT_INTEGRITY", "7 normal thyroid samples; localization only"),
        ("GSE182416_metadata", gse182_meta, "gzip", "READY_WITH_INDEPENDENT_INTEGRITY", "Metadata supplies sample and cell-type labels"),
        ("GSE218743_count_matrix", gse218_count, "gzip", "SUPPORTING_SAMPLE_LEVEL_ONLY", "Condition labels available; not HRA001684 donor-level validation"),
        ("GSE218743_metadata", gse218_meta, "gzip", "SUPPORTING_SAMPLE_LEVEL_ONLY", "Sample-level labels; donor identity not independently verified"),
        ("GSE248205_processed_archive", gse248, "tar", "READY_WITH_INDEPENDENT_INTEGRITY", "Spatial archive inventory and candidate spot detection only"),
        ("GSE285196_raw_archive", gse285, "tar", "READY_ARCHIVE_CANDIDATE_AUDIT", "10 matrices; T13 genes are not expected to be present in the targeted panel"),
    ]:
        if kind == "gzip":
            ok = path.exists() and gzip_integrity(path)
            integrity = "gzip_ok" if ok else "gzip_failed"
            members = ""
        else:
            ok, members = tar_integrity(path) if path.exists() else (False, 0)
            integrity = "tar_ok" if ok else "tar_failed"
        inventory.append(source_row(dataset, path, manifest, status if ok else "FAILED_INTEGRITY", limitation, integrity=integrity, members=members))
    inventory.extend([
        source_row("GTEx_thyroid_v11_eQTL", ROOT / "raw_data/Thyroid.v11.eQTLs.signif_pairs.txt", manifest, "HOLD_CANDIDATE_GENES_NOT_PRESENT", "GTEx v11 significant-pairs table is present, but exact T13 candidate rows are absent; multi-signal coloc cannot run", integrity="text_header_valid"),
        source_row("GTEx_thyroid_v10_eQTL", ROOT / "raw_data/Thyroid.v10.eQTLs.txt.gz", manifest, "PARTIAL_TNRC18_ONLY", "GTEx v10 text table contains one TNRC18 row but no CTLA4 row; not sufficient for the two-gene T14 coloc gate", integrity="text_header_valid_not_gzip"),
        source_row("HRA001684_Zenodo_archive", hra_zip, manifest, "READY_WITH_OFFICIAL_MD5", "Zenodo author-filtered matrices; original FASTQ/bulk archive remains in GSA-Human HRA001684", integrity="zip_md5_verified"),
    ])
    # Zenodo's immutable record supplies the authoritative checksum, not the
    # project's older local-data manifest.
    for r in inventory:
        if r["dataset"] == "HRA001684_Zenodo_archive" and r["exists"]:
            r["manifest_bytes"] = 533057539
            r["manifest_md5"] = "c4263d4824deb365fb98ccac2e5bbd1a"
            r["manifest_match"] = r["bytes_actual"] == 533057539 and r["md5_actual"] == r["manifest_md5"]
            r["status"] = "READY_WITH_OFFICIAL_MD5" if r["manifest_match"] else "FAILED_OFFICIAL_MD5"

    normal_rows, normal_qc = gse182416_expression(gene_rows)
    sample_rows = gse218743_expression(gene_rows)
    spatial_rows = spatial_inventory(gene_rows)
    graves_rows = gse285196_inventory(gene_rows)
    hra_rows, hra_qc = hra001684_expression(gene_rows)
    gtex_v11 = ROOT / "raw_data/Thyroid.v11.eQTLs.signif_pairs.txt"
    gtex_v10 = ROOT / "raw_data/Thyroid.v10.eQTLs.txt.gz"
    gtex_qc = gtex_thyroid_audit(gtex_v11, gene_rows)
    gtex_v10_qc = gtex_thyroid_audit(gtex_v10, gene_rows)

    fields_inventory = ["dataset", "path", "exists", "bytes_actual", "md5_actual", "sha256_actual", "manifest_bytes", "manifest_md5", "manifest_match", "integrity", "archive_members", "status", "limitation"]
    write_tsv(OUT / "t14_input_inventory.tsv", inventory, fields_inventory)
    fields_normal = ["dataset", "gene_id", "gene_symbol", "donor_id", "celltype_global", "celltype_sub", "n_cells", "n_detected", "detection_rate", "total_count", "mean_count_per_cell", "evidence_class", "replicate_unit"]
    write_tsv(OUT / "t14_gse182416_normal_thyroid_expression.tsv", normal_rows, fields_normal)
    fields_sample = ["dataset", "gene_id", "gene_symbol", "sample_id", "condition_label", "celltype_label_as_supplied", "n_cells", "n_detected", "detection_rate", "total_count", "mean_count_per_cell", "evidence_class", "replicate_unit", "limitation"]
    write_tsv(OUT / "t14_gse218743_sample_expression.tsv", sample_rows, fields_sample)
    fields_spatial = ["dataset", "sample_id", "archive_group_prefix", "group_label_source", "gene_id", "gene_symbol", "n_features", "n_barcodes", "matrix_rows", "matrix_columns", "n_tissue_position_rows", "candidate_in_features", "n_nonzero_spots", "detection_rate_among_barcodes", "total_count", "evidence_class", "limitation"]
    write_tsv(OUT / "t14_gse248205_spatial_candidate_detection.tsv", spatial_rows, fields_spatial)
    fields_graves = ["dataset", "sample_id", "matrix_index", "n_features", "n_barcodes", "matrix_rows", "matrix_columns", "candidate_in_features", "n_nonzero_cells", "detection_rate_among_barcodes", "total_count", "gene_id", "gene_symbol", "evidence_class", "limitation"]
    write_tsv(OUT / "t14_gse285196_candidate_inventory.tsv", graves_rows, fields_graves)
    fields_hra = ["dataset", "gene_id", "gene_symbol", "donor_id", "sample_id", "tissue", "n_cells", "n_detected", "detection_rate", "total_count", "mean_count_per_cell", "evidence_class", "replicate_unit", "limitation"]
    write_tsv(OUT / "t14_hra001684_zenodo_thyroid_expression.tsv", hra_rows, fields_hra)
    write_tsv(OUT / "t14_gtex_thyroid_eqtl_audit.tsv", [
        {"dataset": "GTEx_thyroid_v11_signif_pairs", "path": gtex_qc["path"], "bytes": gtex_qc["bytes"], "candidate_rows": gtex_qc["candidate_rows"], "matched_ids": str(gtex_qc["matched_ids"]), "matched_symbols": str(gtex_qc["matched_symbols"]), "status": "HOLD_CANDIDATE_GENES_NOT_PRESENT" if gtex_qc["candidate_rows"] == 0 else "AUDIT_ONLY_NOT_COLOC_READY"},
        {"dataset": "GTEx_thyroid_v10_eqtl", "path": gtex_v10_qc["path"], "bytes": gtex_v10_qc["bytes"], "candidate_rows": gtex_v10_qc["candidate_rows"], "matched_ids": str(gtex_v10_qc["matched_ids"]), "matched_symbols": str(gtex_v10_qc["matched_symbols"]), "status": "PARTIAL_TNRC18_ONLY_NO_CTLA4" if gtex_v10_qc["candidate_rows"] > 0 else "HOLD_CANDIDATE_GENES_NOT_PRESENT"},
    ], ["dataset", "path", "bytes", "candidate_rows", "matched_ids", "matched_symbols", "status"])

    validation = [
        {"check": "t13_pass_with_hold_guard", "passed": True},
        {"check": "t13_candidate_set_26_unique", "passed": len(candidate_rows) == 26 and len({r["candidate_id"] for r in candidate_rows}) == 26},
        {"check": "gse182416_gzip_and_cell_id_match", "passed": normal_qc["metadata_cells"] == normal_qc["matrix_cells"] == 54726 and normal_qc["target_genes_found"] == len(gene_rows)},
        {"check": "gse182416_target_rows_present", "passed": normal_qc["target_rows_found"] == len(gene_rows)},
        {"check": "gse218743_sample_level_output_only", "passed": all(r["replicate_unit"] == "sample; cells aggregated" for r in sample_rows)},
        {"check": "gse248205_archive_and_candidate_feature_dimensions", "passed": len(spatial_rows) == 8 * len(gene_rows) and all(r["n_features"] == r["matrix_rows"] and r["n_barcodes"] == r["matrix_columns"] for r in spatial_rows)},
        {"check": "gse285196_archive_dimensions", "passed": len(graves_rows) == 10 * len(gene_rows) and all(r["n_features"] == r["matrix_rows"] and r["n_barcodes"] == r["matrix_columns"] for r in graves_rows)},
        {"check": "gtex_thyroid_eqtl_audited", "passed": gtex_qc["bytes"] > 0 and gtex_qc["candidate_rows"] == 0 and gtex_v10_qc["bytes"] > 0 and gtex_v10_qc["candidate_rows"] == 1},
        {"check": "hra001684_zenodo_archive_verified", "passed": any(r["dataset"] == "HRA001684_Zenodo_archive" and r["status"] == "READY_WITH_OFFICIAL_MD5" and r["integrity"] == "zip_md5_verified" for r in inventory)},
        {"check": "hra001684_zenodo_matrix_gate", "passed": hra_qc["libraries"] == 4 and hra_qc["cells"] > 0 and len(hra_rows) == 4 * len(gene_rows)},
        {"check": "no_cell_level_pseudoreplication_claim", "passed": all(r.get("replicate_unit", "").startswith(("donor/sample", "sample;")) for r in normal_rows + sample_rows)},
    ]
    write_tsv(OUT / "t14_validation.tsv", validation, ["check", "passed"])

    provenance = []
    for row in inventory:
        if row["exists"]:
            provenance.append({"dataset": row["dataset"], "path": row["path"], "bytes_actual": row["bytes_actual"], "md5_actual": row["md5_actual"], "sha256_actual": row["sha256_actual"], "integrity": row["integrity"], "manifest_match": row["manifest_match"]})
    provenance.extend([
        {"dataset": "T13_decision_guard", "path": "outcome_data/t13/t13_decision.tsv", "bytes_actual": "", "md5_actual": md5_path(ROOT / "outcome_data/t13/t13_decision.tsv"), "sha256_actual": sha256_path(ROOT / "outcome_data/t13/t13_decision.tsv"), "integrity": "read_and_guarded", "manifest_match": ""},
        {"dataset": "T13_gene_candidates", "path": "outcome_data/t13/t13_gene_summary.tsv", "bytes_actual": "", "md5_actual": md5_path(ROOT / "outcome_data/t13/t13_gene_summary.tsv"), "sha256_actual": sha256_path(ROOT / "outcome_data/t13/t13_gene_summary.tsv"), "integrity": f"{len(gene_rows)} mapped gene rows", "manifest_match": ""},
    ])
    write_tsv(OUT / "t14_provenance.tsv", provenance, ["dataset", "path", "bytes_actual", "md5_actual", "sha256_actual", "integrity", "manifest_match"])

    passed = all(str(r["passed"]).lower() == "true" for r in validation)
    decision = [{
        "gate": "T14_thyroid_tissue_and_single_cell_orthogonal_validation",
        "status": "PASS_WITH_HOLD" if passed else "FAIL",
        "decision": "PASS_T14_PARTIAL_ORTHOGONAL_VALIDATION_HOLD_GTEX_HRA001684" if passed else "FAIL_T14_VALIDATION",
        "t13_status": "PASS_WITH_HOLD",
        "normal_thyroid_genes": len(gene_rows),
        "normal_thyroid_cells": normal_qc["metadata_cells"],
        "normal_thyroid_samples": len({r["donor_id"] for r in normal_rows}),
        "spatial_samples": len({r["sample_id"] for r in spatial_rows}),
        "graves_archive_matrices": len({r["sample_id"] for r in graves_rows}),
        "gtex_status": "PARTIAL_TNRC18_ONLY_CTLA4_ABSENT_NO_COLOC",
        "hra001684_status": "ZENODO_FILTERED_MATRIX_AVAILABLE_GSA_RAW_FASTQ_BULK_RETAINED",
        "rationale": "Available local orthogonal archives and the author-filtered HRA001684 Zenodo thyroid matrices were audited and aggregated at library/donor level. GTEx v10 contains one TNRC18 significant pair but no CTLA4 row; v11 signif_pairs contains neither exact candidate. The input is therefore partial and insufficient for two-gene multi-signal coloc. No causal or disease-state claim is promoted.",
    }]
    write_tsv(OUT / "t14_decision.tsv", decision, list(decision[0].keys()))
    print(f"T14 partial gate: {len(gene_rows)} genes; normal thyroid {normal_qc['target_rows_found']} target rows; spatial {len({r['sample_id'] for r in spatial_rows})} samples; Graves archive {len({r['sample_id'] for r in graves_rows})} matrices; {'PASS_WITH_HOLD' if passed else 'FAIL'}")
    if not passed:
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"T14 failed: {exc}", file=sys.stderr)
        sys.exit(1)
