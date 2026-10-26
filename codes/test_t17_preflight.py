"""Small offline fixtures for parser edge cases; does not touch project data."""
import importlib.util
from pathlib import Path
import tempfile
import gzip
import unittest

spec = importlib.util.spec_from_file_location('t17', Path(__file__).with_name('50_preflight_t17_bbj.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def record():
    return dict(CHR='1', POS='751343', SNPID='1_751343_T_A', Allele1='T', Allele2='A',
                AC_Allele2='59762.28125', AF_Allele2='0.140648', N='212453', BETA='0.04118',
                SE='0.04874', **{'p.value': '0.398', 'AF.Cases': '0.1439',
                               'AF.Controls': '0.14046', 'Rsq': '0.81495',
                               'MAC': '59762', 'Is.SPA.converge': '1'})


class QC(unittest.TestCase):
    def test_valid(self):
        self.assertEqual(m.row_errors(record()), [])

    def test_nonfinite_and_invalid(self):
        r = record(); r.update(BETA='nan', SE='0', N='212454', Rsq='0.6')
        self.assertEqual(set(m.row_errors(r)), {'BETA', 'SE', 'N', 'Rsq'})

    def test_p_underflow_retained(self):
        r = record(); r['p.value'] = '0'
        self.assertEqual(m.row_errors(r), [])

    def test_spa_not_wald_p(self):
        r = record(); r['p.value'] = '1e-12'; r['Is.SPA.converge'] = '0'
        self.assertEqual(m.row_errors(r), [])  # flagged for review, not replaced with Wald P

    def test_allele_and_variant_id(self):
        r = record(); r['Allele2'] = 'T'
        self.assertEqual(set(m.row_errors(r)), {'alleles', 'SNPID'})

    def test_variant_specific_n(self):
        r = record(); r['N'] = '212451'
        self.assertEqual(m.row_errors(r), [])

    def test_symbolic_not_malformed(self):
        r = record(); r.update(Allele2='<CN0>', SNPID='1_751343_T_<CN0>')
        self.assertEqual(m.row_errors(r), [])

    def test_duplicate_and_unsorted(self):
        rows = [record(), record(), dict(record(), POS='10', SNPID='1_10_T_A')]
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'GD.auto.rsq07.mac10.txt.gz'
            with gzip.open(p, 'wt') as f:
                f.write(' '.join(rows[0]) + '\n')
                for r in rows:
                    f.write(' '.join(r.values()) + '\n')
            s = m.scan(p)
        self.assertEqual(s['counts']['duplicate_variant_id'], 1)
        self.assertEqual(s['counts']['position_decreases'], 1)


if __name__ == '__main__':
    unittest.main()
