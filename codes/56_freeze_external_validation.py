#!/usr/bin/env python3
"""Freeze protocol/targets before extracting external target associations."""
import csv
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outcome_data/external_validation/plan_20260904'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    assert not OUT.exists(),'Refusing to overwrite frozen plan'
    cfgpath=ROOT/'config/t16_t17_validation_plan.json';cfg=json.loads(cfgpath.read_text())
    targetpath=ROOT/cfg['target_source']
    with targetpath.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==11
    excluded=[r['signal_unit_id'] for r in rows if r['five_class_primary']=='unresolved']
    assert excluded==['GSU058']
    rows=[r for r in rows if r['five_class_primary']!='unresolved']
    assert len(rows)==10 and len({r['representative_variant'] for r in rows})==10
    assert cfg['planned_tests']==20 and cfg['method_review']=='PASS_WITH_EXPLICIT_APPROXIMATION_LIMITS'
    fields=['signal_unit_id','five_class_primary','representative_variant','ref','alt','representative_nearest_genes']
    OUT.mkdir(parents=True)
    (OUT/'plan.json').write_text(cfgpath.read_text())
    with (OUT/'targets.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        w.writerows({k:r[k] for k in fields} for r in rows)
    inputs=[cfgpath,targetpath,ROOT/'config/analysis.yml',ROOT/'config/traits.tsv',
        Path(__file__).resolve(),ROOT/'codes/57_run_external_validation.py',OUT/'plan.json',OUT/'targets.tsv',
        ROOT/'raw_data/external_gwas/bbj_2021_hashimoto/GCST90018635_buildGRCh37.tsv.gz',
        ROOT/'raw_data/external_gwas/bbj_2019_graves/GD.auto.rsq07.mac10.txt.gz',
        *[ROOT/s['discovery_input'] for s in cfg['traits'].values()]]
    result={'frozen_at_utc':datetime.now(timezone.utc).isoformat(),'preregistration':False,
        'scope':'execution-time analytical freeze; prior coverage/discovery results known',
        'sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs}}
    (OUT/'freeze.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
