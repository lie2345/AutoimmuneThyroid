#!/usr/bin/env python3
"""Full offline integrity/QC scan for BBJ Hashimoto GCST90018635.

This program never harmonises alleles, filters variants, runs replication, or
manufactures per-variant sample sizes.  It validates a complete official gzip
source, then writes a new immutable QC run directory only after the scan.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3
import tempfile

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw_data/external_gwas/bbj_2021_hashimoto"
OUT_ROOT = ROOT / "outcome_data/t16/full_qc"
FILENAME = "GCST90018635_buildGRCh37.tsv.gz"
EXPECTED_BYTES = 526_237_307
EXPECTED_MD5 = "91fc15fd74923c4f01ce31f2c87cb2ef"
REQUIRED_COLUMNS = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "effect_allele_frequency", "beta", "standard_error", "p_value", "variant_id",
]
AUTOSOMES = {str(i) for i in range(1, 23)}
ALLOWED_CHROMOSOMES = AUTOSOMES | {"X", "Y", "MT", "M"}
SYMBOLIC = re.compile(r"<[A-Za-z0-9_.:+-]+>")
VARIABLE = re.compile(r"[A-Za-z0-9_.:+-]+")


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def classify_allele(value: str) -> str:
    """Classify source text form; noncanonical form is not biological validation."""
    if not value or any(c.isspace() for c in value):
        return "invalid"
    if set(value.upper()) <= set("ACGT"):
        return "canonical"
    if SYMBOLIC.fullmatch(value):
        return "symbolic"
    if VARIABLE.fullmatch(value):
        return "variable"
    return "invalid"


def valid_number(value: str, predicate) -> bool:
    try:
        return math.isfinite(float(value)) and predicate(float(value))
    except (TypeError, ValueError):
        return False


def row_errors(row: dict[str, str]) -> list[str]:
    errors = []
    if not row["base_pair_location"].isdigit() or int(row["base_pair_location"]) <= 0:
        errors.append("base_pair_location")
    if not valid_number(row["effect_allele_frequency"], lambda x: 0 <= x <= 1):
        errors.append("effect_allele_frequency")
    if not valid_number(row["beta"], lambda x: True):
        errors.append("beta")
    if not valid_number(row["standard_error"], lambda x: x > 0):
        errors.append("standard_error")
    # p=0 is retained as an underflow flag, not rejected or replaced.
    if not valid_number(row["p_value"], lambda x: 0 <= x <= 1):
        errors.append("p_value")
    if "invalid" in (classify_allele(row["effect_allele"]),
                     classify_allele(row["other_allele"])):
        errors.append("allele_encoding")
    if row["effect_allele"].upper() == row["other_allele"].upper():
        errors.append("same_alleles")
    return errors


def is_missing(value: str) -> bool:
    return value.strip().upper() in {"", "NA", "N/A", ".", "NULL"}


def make_database(path: Path) -> sqlite3.Connection:
    db = sqlite3.connect(path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("PRAGMA temp_store=FILE")
    # The primary key is the exact coordinate-plus-allele key.  SQLite keeps
    # this bounded on disk; all records are still parsed from the gzip stream.
    db.execute("""CREATE TABLE exact_key (
        chromosome TEXT NOT NULL, position INTEGER NOT NULL,
        effect_allele TEXT NOT NULL, other_allele TEXT NOT NULL,
        first_row INTEGER NOT NULL, occurrences INTEGER NOT NULL,
        PRIMARY KEY (chromosome, position, effect_allele, other_allele)
    ) WITHOUT ROWID""")
    return db


def scan(path: Path, db_path: Path, progress_every: int = 2_000_000) -> dict:
    """Read every uncompressed row. gzip reaching EOF validates the trailer CRC."""
    counts = collections.Counter()
    chromosomes = collections.Counter()
    errors = collections.Counter()
    examples = []
    first_noncanonical = []
    previous_position: dict[str, int] = {}
    db = make_database(db_path)
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames or []
            if header != REQUIRED_COLUMNS:
                raise ValueError(f"Unexpected schema: {header!r}")
            # Missing N is a downstream limitation of this source schema, not
            # a numeric error to be filled with the Catalog total sample size.
            counts["n_column_absent"] = int("n" not in header and "N" not in header)
            batch = []
            for row_number, row in enumerate(reader, start=2):
                counts["rows"] += 1
                if None in row or any(v is None for v in row.values()):
                    errors["field_count"] += 1
                    if len(examples) < 10:
                        examples.append({"row": row_number, "errors": ["field_count"]})
                    continue
                row_error = row_errors(row)
                for error in row_error:
                    errors[error] += 1
                if row_error and len(examples) < 10:
                    examples.append({"row": row_number, "variant_id": row["variant_id"],
                                     "errors": row_error})

                chrom = row["chromosome"]
                chromosomes[chrom] += 1
                if chrom not in ALLOWED_CHROMOSOMES:
                    counts["unexpected_chromosome_rows"] += 1
                try:
                    pos = int(row["base_pair_location"])
                    if previous_position.get(chrom, pos) > pos:
                        counts["position_decreases_within_chromosome"] += 1
                    previous_position[chrom] = pos
                except ValueError:
                    pos = -1

                ea_type, oa_type = classify_allele(row["effect_allele"]), classify_allele(row["other_allele"])
                if ea_type == oa_type == "canonical":
                    counts["canonical_allele_rows"] += 1
                    if len(row["effect_allele"]) == len(row["other_allele"]) == 1:
                        counts["snv_rows"] += 1
                        if row["effect_allele"].upper() + row["other_allele"].upper() in {"AT", "TA", "CG", "GC"}:
                            counts["palindromic_snv_rows"] += 1
                    else:
                        counts["canonical_indel_or_complex_rows"] += 1
                else:
                    if "symbolic" in (ea_type, oa_type):
                        counts["symbolic_allele_rows"] += 1
                    elif "variable" in (ea_type, oa_type):
                        counts["variable_allele_rows"] += 1
                    if len(first_noncanonical) < 10:
                        first_noncanonical.append({"row": row_number, "effect_allele": row["effect_allele"],
                                                   "other_allele": row["other_allele"], "classes": [ea_type, oa_type]})
                if is_missing(row["variant_id"]):
                    counts["missing_variant_id_rows"] += 1
                if valid_number(row["p_value"], lambda x: x == 0):
                    counts["p_value_zero_rows"] += 1

                # Invalid coordinate, encoding, or same-allele rows remain in
                # the structural error count but are not inserted into an
                # artificial duplicate key.
                if pos > 0 and not {"allele_encoding", "same_alleles"}.intersection(row_error):
                    batch.append((chrom, pos, row["effect_allele"], row["other_allele"], row_number))
                if len(batch) >= 50_000:
                    _flush_keys(db, batch)
                    batch.clear()
                if progress_every and counts["rows"] % progress_every == 0:
                    db.commit()
                    print(f"{path.name}: {counts['rows']:,} rows scanned", flush=True)
            if batch:
                _flush_keys(db, batch)
            db.commit()
            # The context manager has now consumed EOF, so gzip CRC/ISIZE
            # validation has occurred before results are returned.
        duplicate_extra = db.execute("SELECT COALESCE(SUM(occurrences - 1), 0) FROM exact_key WHERE occurrences > 1").fetchone()[0]
        duplicate_distinct = db.execute("SELECT COUNT(*) FROM exact_key WHERE occurrences > 1").fetchone()[0]
        duplicate_examples = [dict(zip(["chromosome", "position", "effect_allele", "other_allele", "occurrences"], r))
                              for r in db.execute("SELECT chromosome, position, effect_allele, other_allele, occurrences FROM exact_key WHERE occurrences > 1 ORDER BY occurrences DESC, chromosome, position LIMIT 10")]
        multi_positions = db.execute("SELECT COUNT(*) FROM (SELECT chromosome, position FROM exact_key GROUP BY chromosome, position HAVING COUNT(*) > 1)").fetchone()[0]
        multi_examples = [dict(zip(["chromosome", "position", "allele_pairs"], r))
                          for r in db.execute("SELECT chromosome, position, COUNT(*) FROM exact_key GROUP BY chromosome, position HAVING COUNT(*) > 1 ORDER BY chromosome, position LIMIT 10")]
    finally:
        db.close()
    counts["duplicate_exact_key_extra_records"] = duplicate_extra
    counts["duplicate_exact_key_distinct_keys"] = duplicate_distinct
    counts["multi_allelic_positions"] = multi_positions
    return {"file": path.name, "header": REQUIRED_COLUMNS, "counts": dict(counts),
            "chromosomes": dict(chromosomes), "errors": dict(errors), "examples": examples,
            "noncanonical_examples": first_noncanonical, "duplicate_examples": duplicate_examples,
            "multi_allelic_position_examples": multi_examples}


def _flush_keys(db: sqlite3.Connection, batch: list[tuple]) -> None:
    db.executemany("""INSERT INTO exact_key (chromosome, position, effect_allele, other_allele, first_row, occurrences)
                      VALUES (?, ?, ?, ?, ?, 1)
                      ON CONFLICT(chromosome, position, effect_allele, other_allele)
                      DO UPDATE SET occurrences = occurrences + 1""", batch)


def promote_input(raw_dir: Path) -> tuple[Path, bool]:
    """Validate published size/MD5 before a .part can become the final name."""
    final = raw_dir / FILENAME
    part = raw_dir / f"{FILENAME}.part"
    candidate = final if final.exists() else part
    if not candidate.exists():
        raise SystemExit("HOLD: no complete source or resumable .part file found")
    if candidate.stat().st_size != EXPECTED_BYTES:
        raise SystemExit(f"HOLD: expected {EXPECTED_BYTES} bytes, found {candidate.stat().st_size} in {candidate.name}")
    observed_md5 = digest(candidate, "md5")
    if observed_md5 != EXPECTED_MD5:
        raise SystemExit(f"HOLD: published MD5 mismatch for {candidate.name}: {observed_md5}")
    promoted = False
    if candidate == part:
        part.replace(final)
        promoted = True
    return final, promoted


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def create_output_dir(requested: Path | None) -> Path:
    if requested is not None:
        if requested.exists():
            raise SystemExit(f"Refusing to overwrite existing QC output: {requested}")
        return requested
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    for suffix in range(1000):
        candidate = OUT_ROOT / (stamp if suffix == 0 else f"{stamp}_{suffix:03d}")
        if not candidate.exists():
            return candidate
    raise SystemExit("Unable to allocate a new QC output directory")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, help="Override default raw source path (must be final .gz or .part)")
    parser.add_argument("--out-dir", type=Path, help="New output directory; refuses an existing path")
    args = parser.parse_args()
    if args.input:
        # Fixtures may use an alternate path, but production promotion is only
        # permitted for the exact registered raw-source location.
        source = args.input
        if not source.exists():
            raise SystemExit(f"HOLD: input does not exist: {source}")
        if source.stat().st_size != EXPECTED_BYTES or digest(source, "md5") != EXPECTED_MD5:
            raise SystemExit("HOLD: override input does not match published source size/MD5")
        promoted = False
    else:
        source, promoted = promote_input(RAW)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".t16_qc_", dir=OUT_ROOT) as staging:
        result = scan(source, Path(staging) / "exact_keys.sqlite")
    observed = set(result["chromosomes"])
    missing_autosomes = sorted(AUTOSOMES - observed, key=int)
    unexpected = sorted(observed - ALLOWED_CHROMOSOMES)
    structural_checks = [
        {"check": "published_file_size", "status": "PASS", "detail": str(EXPECTED_BYTES)},
        {"check": "published_md5", "status": "PASS", "detail": EXPECTED_MD5},
        {"check": "gzip_full_stream_crc", "status": "PASS", "detail": f"all {result['counts']['rows']} rows reached EOF"},
        {"check": "exact_expected_schema", "status": "PASS", "detail": "tab-delimited 9-column pre-GWAS-SSF header"},
        {"check": "all_rows_have_required_fields", "status": "PASS" if not result["errors"].get("field_count") else "FAIL", "detail": str(result["errors"].get("field_count", 0))},
        {"check": "allowed_chromosome_labels", "status": "PASS" if not unexpected else "FAIL", "detail": json.dumps({"unexpected": unexpected})},
        {"check": "all_autosomes_observed", "status": "PASS" if not missing_autosomes else "FAIL", "detail": json.dumps({"missing": missing_autosomes})},
        {"check": "finite_valid_numeric_fields", "status": "PASS" if not any(result["errors"].get(k) for k in ["effect_allele_frequency", "beta", "standard_error", "p_value", "base_pair_location"]) else "FAIL", "detail": json.dumps(result["errors"])},
        {"check": "valid_allele_encodings", "status": "PASS" if not result["errors"].get("allele_encoding") else "FAIL", "detail": str(result["errors"].get("allele_encoding", 0))},
        {"check": "effect_and_other_alleles_differ", "status": "PASS" if not result["errors"].get("same_alleles") else "FAIL", "detail": str(result["errors"].get("same_alleles", 0))},
    ]
    structural_complete = all(c["status"] == "PASS" for c in structural_checks)
    holds = [
        "No per-variant N column is present; Catalog 537 cases/172656 controls is cohort-level metadata and was not filled into rows.",
        "variant_id missingness is retained; coordinate-plus-allele keys were used only for QC, not to create cleaned IDs.",
        "Duplicate coordinate-plus-allele keys and multi-allelic positions require an explicit downstream selection/normalisation rule.",
        "Symbolic and variable allele text forms (including forms such as NA or -) are format flags only; they are not evidence of biological validity and require a downstream-compatible harmonisation path.",
        "No liftover, effect alignment, ancestry-matched LD, replication test, or European independent AIHT source was executed.",
    ]
    decision = {
        "completed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "stage": "T16 BBJ Hashimoto full acquisition/QC only",
        "status": "PASS_WITH_HOLD" if structural_complete else "HOLD_QC_REVIEW",
        "structural_complete": structural_complete,
        "input_promoted_from_part_after_size_md5": promoted,
        "rows_scanned": result["counts"]["rows"],
        "expected_catalog_cohort": {"cases": 537, "controls": 172656, "total": 173193,
                                    "note": "Repeated YAML sample entries were not summed."},
        "genome_build": "GRCh37", "is_harmonised": False, "is_sorted": False,
        "replication_executed": False, "liftover_executed": False,
        "holds": holds,
    }
    output = create_output_dir(args.out_dir)
    output.mkdir(parents=True)
    (output / "t16_bbj_file_qc.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_tsv(output / "t16_bbj_validation.tsv", structural_checks)
    (output / "t16_bbj_decision.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    inputs = [source, Path(__file__).resolve(), ROOT / "outcome_data/t16/sources/GCST90018635.json",
              ROOT / "outcome_data/t16/sources/GCST90018635_meta.yaml"]
    provenance = []
    for item in inputs:
        provenance.append({"path": str(item.relative_to(ROOT)), "bytes": item.stat().st_size,
                           "sha256": digest(item, "sha256"), "md5": ""})
    provenance.append({"path": "reference:published_file", "bytes": EXPECTED_BYTES,
                       "sha256": "not_published", "md5": EXPECTED_MD5})
    write_tsv(output / "t16_bbj_provenance.tsv", provenance)
    print(json.dumps({"status": decision["status"], "rows_scanned": decision["rows_scanned"],
                      "output_dir": str(output.relative_to(ROOT))}, indent=2), flush=True)
    if not structural_complete:
        raise SystemExit("HOLD: full scan completed but structural validation requires review")


if __name__ == "__main__":
    main()
