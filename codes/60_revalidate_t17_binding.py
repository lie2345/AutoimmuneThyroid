#!/usr/bin/env python3
"""Re-run unchanged T17 independent audit with before/after input hashes."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'outcome_data/t17/harmonization_resume_20260904'
OUT=ROOT/'outcome_data/t17/revalidation_bound_20260904'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    assert not OUT.exists();OUT.mkdir()
    names=['harmonization.sqlite','bbj_finngen_graves.aligned.grch38.tsv.gz','decision.json']
    before={str((SOURCE/n).relative_to(ROOT)):sha(SOURCE/n) for n in names}
    for n in names:(OUT/n).symlink_to(SOURCE/n)
    subprocess.run([sys.executable,str(ROOT/'codes/54_validate_t17_harmonization.py'),'--output',str(OUT.relative_to(ROOT))],check=True)
    audit=json.loads((OUT/'independent_validation.json').read_text());assert audit['status']=='PASS'
    for n in names:assert sha(SOURCE/n)==before[str((SOURCE/n).relative_to(ROOT))]
    (OUT/'verified_binding.json').write_text(json.dumps({'status':'PASS','verified_input_sha256':before,
        'independent_audit_sha256':sha(OUT/'independent_validation.json')},indent=2)+'\n')
    print('T17_BOUND_FULL_AUDIT_PASS',flush=True)


if __name__=='__main__':main()
