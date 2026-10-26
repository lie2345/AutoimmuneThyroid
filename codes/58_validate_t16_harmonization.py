#!/usr/bin/env python3
"""Independent raw-to-table and full exported-row audit for T16."""
import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='outcome_data/t16/harmonization_20260904')
    args=ap.parse_args();out=(ROOT/args.output).resolve()
    assert out.is_relative_to(ROOT/'outcome_data/t16') and not (out/'independent_validation.json').exists()
    decision=json.loads((out/'decision.json').read_text());assert decision['status'].startswith('PASS')
    inputs=[out/'harmonization.sqlite',out/'bbj_finngen_aiht.aligned.grch38.tsv.gz',out/'decision.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    db=sqlite3.connect(f'file:{out}/harmonization.sqlite?mode=ro',uri=True)
    cur=iter(db.execute('SELECT id,c37,p37,ea37,oa37,beta,se,p,eaf FROM v ORDER BY id'))
    target=next(cur,None);counts=collections.Counter();errors=collections.Counter()
    with gzip.open(ROOT/'raw_data/external_gwas/bbj_2021_hashimoto/GCST90018635_buildGRCh37.tsv.gz','rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        for i,r in enumerate(reader,1):
            counts['raw_rows_scanned']+=1
            if target is None or i!=target[0]:continue
            observed=(i,int(r['chromosome']),int(r['base_pair_location']),r['effect_allele'],r['other_allele'],
                float(r['beta']),float(r['standard_error']),r['p_value'],float(r['effect_allele_frequency']))
            if observed!=target:errors['retained_raw_mismatch']+=1
            counts['retained_raw_verified']+=1;target=next(cur,None)
    if target is not None:errors['unmatched_retained_rows']+=1
    print('T16 retained rows verified against source',flush=True)
    sql="""SELECT v.id,v.c38,v.p38,v.p37,v.ea37,v.oa37,v.beta,v.se,v.p,v.eaf,v.strand,v.status,
        f.vid,f.ref,f.alt,f.beta,f.se,f.pval,f.eaf FROM v JOIN fg f ON f.c=v.c38 AND f.p=v.p38
        WHERE v.status IN ('ALIGNED_DIRECT','ALIGNED_SWAP') ORDER BY v.c38,v.p38"""
    comp={'A':'T','T':'A','C':'G','G':'C'};previous=(0,0)
    with gzip.open(out/'bbj_finngen_aiht.aligned.grch38.tsv.gz','rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        for r,x in itertools.zip_longest(reader,db.execute(sql)):
            if r is None or x is None:errors['output_row_count']+=1;continue
            i,c,p,p37,ea,oa,b,se,sp,eaf,strand,status,vid,fr,fa,fb,fs,fp,fe=x
            source_ea,source_oa=ea,oa
            if strand=='-':ea,oa=comp[ea],comp[oa]
            direct=(ea,oa)==(fa,fr);swap=(ea,oa)==(fr,fa)
            if direct==swap:errors['allele_orientation']+=1
            expected_status='ALIGNED_SWAP' if swap else 'ALIGNED_DIRECT'
            if expected_status!=status:errors['sql_status_mismatch']+=1
            numeric={'source_row':i,'chromosome':c,'position_grch38':p,'bbj_beta_aligned':-b if swap else b,
                'bbj_se':se,'bbj_eaf_aligned':1-eaf if swap else eaf,'finngen_beta':fb,'finngen_se':fs,
                'finngen_eaf':fe,'source_position_grch37':p37,'bbj_beta_source':b,'bbj_eaf_source':eaf}
            strings={'discovery_variant_id':vid,'effect_allele':fa,'other_allele':fr,'bbj_p':sp,'finngen_p':fp,
                'source_effect_allele':source_ea,'source_other_allele':source_oa,'lift_strand':strand,'alignment':status,'bbj_n_total':'NA'}
            for k,v in numeric.items():
                if float(r[k])!=v:errors['export_mismatch:'+k]+=1
            for k,v in strings.items():
                if r[k]!=v:errors['export_mismatch:'+k]+=1
            if (c,p)<=previous:errors['target_sort_or_unique']+=1
            if source_ea+source_oa in ['AT','TA','CG','GC']:errors['palindrome']+=1
            if c==6 and (25000000<=p<=34000000 or 25000000<=p37<=34000000):errors['hla']+=1
            previous=(c,p);counts['aligned_rows_verified']+=1;counts[status+':'+strand]+=1
    if counts['raw_rows_scanned']!=13428945:errors['raw_count']+=1
    if counts['retained_raw_verified']!=decision['clean_records']:errors['retained_count']+=1
    if counts['aligned_rows_verified']!=decision['aligned_records']:errors['aligned_count']+=1
    bad_roundtrip=db.execute("""SELECT COUNT(*) FROM v LEFT JOIN rv_one r ON r.id=v.id WHERE v.status LIKE 'ALIGNED_%'
        AND (r.id IS NULL OR r.c!=CAST(v.c37 AS TEXT) OR r.p!=v.p37 OR r.strand!='+' OR r.width!=1)""").fetchone()[0]
    if bad_roundtrip:errors['roundtrip']+=bad_roundtrip
    for p in inputs:
        if sha(p)!=before[str(p.relative_to(ROOT))]:errors['audit_input_changed']+=1
    result={'status':'PASS' if not errors else 'FAIL','counts':dict(counts),'errors':dict(errors),
        'verified_input_sha256':before,
        'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'All retained raw rows, every aligned exported row and exact roundtrip; no statistical tests'}
    (out/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    db.close();assert not errors


if __name__=='__main__':main()
