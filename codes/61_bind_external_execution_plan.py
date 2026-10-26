#!/usr/bin/env python3
"""Execution-only amendment: bind validated files, keep frozen statistics intact."""
import copy
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outcome_data/external_validation/plan_20260904'
OUT=ROOT/'outcome_data/external_validation/execution_plan_20260904'
T16=ROOT/'outcome_data/t16/harmonization_verified_20260904'
T17=ROOT/'outcome_data/t17/revalidation_bound_20260904'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    assert not OUT.exists()
    frozen=json.loads((BASE/'freeze.json').read_text())
    for path,value in frozen['sha256'].items():assert sha(ROOT/path)==value,path
    original=json.loads((BASE/'plan.json').read_text());cfg=copy.deepcopy(original)
    t16_audit=json.loads((T16/'independent_validation.json').read_text())
    t17_binding=json.loads((T17/'verified_binding.json').read_text())
    assert t16_audit['status']==t17_binding['status']=='PASS'
    assert sha(T17/'independent_validation.json')==t17_binding['independent_audit_sha256']
    bindings={**t16_audit['verified_input_sha256'],**t17_binding['verified_input_sha256']}
    for path,value in bindings.items():assert sha(ROOT/path)==value,path
    for key,name in [('aligned_input','bbj_finngen_aiht.aligned.grch38.tsv.gz'),('decision','decision.json'),('independent_audit','independent_validation.json')]:
        cfg['traits']['AIHT'][key]=str((T16/name).relative_to(ROOT))
    cfg['traits']['Graves']['independent_audit']=str((T17/'independent_validation.json').relative_to(ROOT))
    a,b=copy.deepcopy(original),copy.deepcopy(cfg)
    for c in [a,b]:
        for trait in c['traits'].values():
            for key in ['aligned_input','decision','independent_audit']:trait.pop(key)
    assert a==b,'Statistical parameters changed'
    OUT.mkdir()
    (OUT/'plan.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (OUT/'targets.tsv').write_bytes((BASE/'targets.tsv').read_bytes())
    amendment={'created_at_utc':datetime.now(timezone.utc).isoformat(),'type':'execution_only_input_binding_and_routing',
        'reason':'T16 interrupted runs excluded after status-transition correction and process isolation; bind successful fresh output and exact re-audited files before and after extraction',
        'original_plan_sha256':sha(BASE/'plan.json'),'original_freeze_sha256':sha(BASE/'freeze.json'),
        'statistics_unchanged':True,'targets_identical':sha(OUT/'targets.tsv')==sha(BASE/'targets.tsv'),
        'runner_unchanged':True,'verified_input_sha256':bindings,
        't16_final_run':str(T16.relative_to(ROOT)),'t17_fresh_audit':str(T17.relative_to(ROOT))}
    (OUT/'amendment.json').write_text(json.dumps(amendment,indent=2)+'\n')
    entries={**frozen['sha256'],**bindings}
    extra=[OUT/'plan.json',OUT/'targets.tsv',OUT/'amendment.json',BASE/'freeze.json',Path(__file__).resolve(),
        T16/'independent_validation.json',T17/'independent_validation.json',T17/'verified_binding.json']
    for p in extra:entries[str(p.relative_to(ROOT))]=sha(p)
    (OUT/'freeze.json').write_text(json.dumps({'bound_at_utc':datetime.now(timezone.utc).isoformat(),
        'original_statistical_freeze_at_utc':frozen['frozen_at_utc'],'sha256':entries},indent=2)+'\n')
    print(json.dumps(amendment,indent=2))


if __name__=='__main__':main()
