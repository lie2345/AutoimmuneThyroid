"""Offline unit tests for T16 source-format QC; does not touch project data."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("t16", Path(__file__).with_name("52_qc_t16_bbj.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def record(**overrides):
    base = dict(chromosome="1", base_pair_location="751343", effect_allele="A", other_allele="T",
                effect_allele_frequency="0.1477", beta="-0.0374", standard_error="0.0901",
                p_value="0.6783", variant_id="NA")
    base.update(overrides)
    return base


class T16QC(unittest.TestCase):
    def test_legal_allele_classes(self):
        self.assertEqual(m.classify_allele("AC"), "canonical")
        self.assertEqual(m.classify_allele("<DEL>"), "symbolic")
        self.assertEqual(m.classify_allele("D"), "variable")
        self.assertEqual(m.classify_allele("NA"), "variable")  # source text form, not biological validation
        self.assertEqual(m.classify_allele("-"), "variable")
        self.assertEqual(m.classify_allele("A T"), "invalid")

    def test_numeric_validation_and_zero_p(self):
        self.assertEqual(m.row_errors(record(p_value="0")), [])
        bad = record(beta="nan", standard_error="0", effect_allele_frequency="1.1", p_value="2")
        self.assertEqual(set(m.row_errors(bad)), {"beta", "standard_error", "effect_allele_frequency", "p_value"})

    def test_scan_flags_duplicates_multiallelic_and_missing_ids(self):
        rows = [record(), record(), record(effect_allele="C", other_allele="T"),
                record(chromosome="2", base_pair_location="10", effect_allele="<DEL>", other_allele="A", variant_id="rs1")]
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "fixture.tsv.gz"
            with gzip.open(source, "wt", newline="") as handle:
                handle.write("\t".join(m.REQUIRED_COLUMNS) + "\n")
                for row in rows:
                    handle.write("\t".join(row[k] for k in m.REQUIRED_COLUMNS) + "\n")
            got = m.scan(source, Path(tmp) / "keys.sqlite", progress_every=0)
        self.assertEqual(got["counts"]["rows"], 4)
        self.assertEqual(got["counts"]["missing_variant_id_rows"], 3)
        self.assertEqual(got["counts"]["duplicate_exact_key_extra_records"], 1)
        self.assertEqual(got["counts"]["multi_allelic_positions"], 1)
        self.assertEqual(got["counts"]["symbolic_allele_rows"], 1)
        self.assertEqual(got["counts"]["n_column_absent"], 1)

    def test_missing_id_is_flag_not_numeric_or_allele_error(self):
        self.assertEqual(m.row_errors(record(variant_id=".")), [])

    def test_same_effect_and_other_allele_is_structural_error(self):
        self.assertEqual(m.row_errors(record(effect_allele="A", other_allele="a")), ["same_alleles"])


if __name__ == "__main__":
    unittest.main()
