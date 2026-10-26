import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('t16',HERE/'55_harmonize_t16_bbj.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
CFG=json.loads((HERE.parent/'config/t16_harmonization.json').read_text())


class SQLAlignment(unittest.TestCase):
    def test_all_strand_swap_combinations(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp);source=d/'fg.gz'
            with gzip.open(source,'wt') as f:
                f.write('chrom\tpos\tref\talt\tbeta\tse\tp\teaf\tvariant_id\n')
                for i,(ref,alt) in enumerate([('C','A'),('A','C'),('G','T'),('T','G')],1):
                    f.write(f'1\t{100+i}\t{ref}\t{alt}\t.3\t.1\t.001\t.2\tchr1_{100+i}_{ref}_{alt}\n')
            db=sqlite3.connect(':memory:')
            db.execute('CREATE TABLE v(id INTEGER,c37 INTEGER,p37 INTEGER,ea37 TEXT,oa37 TEXT,beta REAL,se REAL,p TEXT,eaf REAL,status TEXT,c38 INTEGER,p38 INTEGER,strand TEXT,ea38 TEXT,oa38 TEXT)')
            for i in range(1,5):
                strand='+' if i<=2 else '-';ea,oa=('A','C') if i<=2 else ('T','G')
                db.execute('INSERT INTO v VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(i,1,10+i,'A','C',.2,.1,'0.1234',.1,'READY',1,100+i,strand,ea,oa))
            m.align(db,d,dict(CFG,discovery_input=str(source)))
            with gzip.open(d/'bbj_finngen_aiht.aligned.grch38.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),4)
            self.assertEqual([r['alignment'] for r in rows],['ALIGNED_DIRECT','ALIGNED_SWAP']*2)
            self.assertEqual([float(r['bbj_beta_aligned']) for r in rows],[.2,-.2,.2,-.2])
            self.assertEqual([float(r['bbj_eaf_aligned']) for r in rows],[.1,.9,.1,.9])
            self.assertTrue(all(r['bbj_n_total']=='NA' and r['bbj_p']=='0.1234' for r in rows))
            db.close()


if __name__=='__main__':unittest.main()
