import importlib.util
import json
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('t16', HERE/'55_harmonize_t16_bbj.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
CFG = json.loads((HERE.parent/'config/t16_harmonization.json').read_text())


class T16HarmonizationTests(unittest.TestCase):
    def test_effect_other_semantics_and_swap(self):
        self.assertEqual(m.primary_reason(1, 100, 'A', 'C', .005, 0, .1, .5, CFG, 'source'), None)
        self.assertEqual(m.primary_reason(1, 100, 'A', 'T', .1, 0, .1, .5, CFG, 'source'), 'PALINDROMIC_SNV')
        # EA/OA on a negative lift strand become T/G: direct to FG ALT/REF.
        ea, oa = 'A'.translate(m.COMP), 'C'.translate(m.COMP)
        self.assertEqual((ea, oa), ('T', 'G'))
        self.assertEqual(('T', 'G'), ('T', 'G'))
        # If transformed EA/OA equals FG REF/ALT, alignment is swap and beta/EAF flip only then.
        self.assertEqual((-(-.2), 1-.1), (.2, .9))

    def test_load_bed_multiple_mappings_and_actual_negative_roundtrip(self):
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            (d/'forward.chain').write_text('chain 1000 chr1 1000 + 100 200 chr1 1000 - 300 400 1\n100\n\n')
            (d/'reverse.chain').write_text('chain 1000 chr1 1000 + 600 700 chr1 1000 - 800 900 1\n100\n\n')
            (d/'source.bed').write_text('chr1\t109\t110\t1\t0\t+\n')
            exe = str(HERE.parent/'tools/liftOver')
            subprocess.run([exe,'-minMatch=1.0','-multiple','-noSerial',str(d/'source.bed'),str(d/'forward.chain'),str(d/'mapped.bed'),str(d/'mapped.unmapped')], check=True, capture_output=True)
            subprocess.run([exe,'-minMatch=1.0','-multiple','-noSerial',str(d/'mapped.bed'),str(d/'reverse.chain'),str(d/'returned.bed'),str(d/'returned.unmapped')], check=True, capture_output=True)
            self.assertEqual((d/'mapped.bed').read_text(), 'chr1\t690\t691\t1\t0\t-\n')
            self.assertEqual((d/'returned.bed').read_text(), (d/'source.bed').read_text())
            db = sqlite3.connect(':memory:')
            (d/'many.bed').write_text('chr1\t99\t100\t1\t0\t+\nchr1\t99\t100\t1\t0\t-\n')
            m.load_bed(db, 'fw', d/'many.bed')
            self.assertEqual(db.execute('SELECT COUNT(*) FROM fw_one').fetchone()[0], 0)

    def test_lift_status_end_to_end_reaches_ready(self):
        """Regression for status assignment after unmapped defaults are set."""
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            (d/'forward.chain').write_text('chain 1000 chr1 1000 + 100 200 chr1 1000 - 300 400 1\n100\n\n')
            (d/'reverse.chain').write_text('chain 1000 chr1 1000 + 600 700 chr1 1000 - 800 900 1\n100\n\n')
            (d/'bbj_clean.grch37.bed').write_text('chr1\t109\t110\t1\t0\t+\n')
            db = sqlite3.connect(':memory:')
            db.execute('''CREATE TABLE v(id INTEGER PRIMARY KEY,c37 INTEGER,p37 INTEGER,ea37 TEXT,oa37 TEXT,beta REAL,se REAL,p TEXT,eaf REAL,status TEXT,c38 INTEGER,p38 INTEGER,strand TEXT,ea38 TEXT,oa38 TEXT)''')
            db.execute("INSERT INTO v(id,c37,p37,ea37,oa37,beta,se,p,eaf,status) VALUES(1,1,110,'A','C',.2,.1,'.5',.1,'CANDIDATE')")
            cfg = dict(CFG, forward_chain=str(d/'forward.chain'), reverse_chain=str(d/'reverse.chain'), liftover_executable=str(HERE.parent/'tools/liftOver'))
            m.lift(db, d, cfg)
            self.assertEqual(db.execute('SELECT status,c38,p38,strand,ea38,oa38 FROM v').fetchone(), ('READY',1,691,'-','T','G'))


if __name__ == '__main__': unittest.main()
