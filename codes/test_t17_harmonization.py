import importlib.util
import json
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t17', HERE/'53_clean_harmonize_t17_bbj.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
CFG = json.loads((HERE.parent/'config/t17_harmonization.json').read_text())


class HarmonizationTests(unittest.TestCase):
    def test_effect_orientation(self):
        for args, expected in [
            (('A','C','+','A','C'),1), (('A','C','+','C','A'),-1),
            (('A','C','-','T','G'),1), (('A','C','-','G','T'),-1),
            (('A','C','+','T','G'),None), (('A','T','+','A','T'),None),
            (('C','G','-','G','C'),None), (('A','AC','+','A','AC'),None),
            (('A','A','+','A','A'),None), (('A','C','?','A','C'),None)]:
            with self.subTest(args=args): self.assertEqual(m.allele_sign(*args),expected)

    def test_filter_boundaries(self):
        def reason(c=1,p=100,a='A',b='C',eaf=.1,beta=0,se=.1,pval=.5):
            return m.primary_reason(c,p,a,b,eaf,beta,se,pval,CFG)
        self.assertIsNone(reason(eaf=.005))
        self.assertEqual(reason(eaf=.00499),'MAF_BELOW_0_005')
        for p in [25000000,34000000]: self.assertEqual(reason(c=6,p=p),'HLA_SOURCE_BUILD')
        for p in [24999999,34000001]: self.assertIsNone(reason(c=6,p=p))
        self.assertEqual(reason(beta=float('nan')),'INVALID_NUMERIC')
        self.assertEqual(reason(se=0),'INVALID_NUMERIC')
        self.assertEqual(reason(a='A',b='T'),'PALINDROMIC_SNV')
        self.assertEqual(reason(b='<CN0>'),'NON_CANONICAL_SNV')
        self.assertEqual(reason(c=23),'NON_AUTOSOME_OR_INVALID_POSITION')

    def test_multiple_mappings_never_choose_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed=Path(tmp)/'mapped.bed'
            bed.write_text('chr1\t99\t100\t1\t0\t+\nchr1\t199\t200\t2\t0\t-\nchr1_alt\t99\t100\t2\t0\t+\n')
            db=sqlite3.connect(':memory:')
            m.load_bed(db,'fw',bed)
            self.assertEqual(db.execute('SELECT * FROM fw_one').fetchall(),[(1,'1',100,'+',1)])
            db.close()

    def test_actual_liftover_negative_strand_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/'forward.chain').write_text('chain 1000 chr1 1000 + 100 200 chr1 1000 - 300 400 1\n100\n\n')
            (d/'reverse.chain').write_text('chain 1000 chr1 1000 + 600 700 chr1 1000 - 800 900 1\n100\n\n')
            (d/'source.bed').write_text('chr1\t109\t110\t1\t0\t+\n')
            exe=str(HERE.parent/'tools/liftOver')
            for source,chain,target in [('source.bed','forward.chain','mapped.bed'),('mapped.bed','reverse.chain','returned.bed')]:
                subprocess.run([exe,'-minMatch=1.0','-multiple','-noSerial',str(d/source),str(d/chain),str(d/target),str(d/(target+'.unmapped'))],check=True,capture_output=True)
            self.assertEqual((d/'mapped.bed').read_text(),'chr1\t690\t691\t1\t0\t-\n')
            self.assertEqual((d/'returned.bed').read_text(),(d/'source.bed').read_text())


if __name__=='__main__': unittest.main()
