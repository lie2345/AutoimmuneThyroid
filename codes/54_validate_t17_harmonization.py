#!/usr/bin/env python3
"""Independent full-stream verification against raw BBJ and aligned output."""
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


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    out=(ROOT/args.output).resolve()
    assert out.is_relative_to(ROOT/'outcome_data/t17')
    assert not (out/'independent_validation.json').exists(), 'Refusing to overwrite audit'
    decision=json.loads((out/'decision.json').read_text())
    assert decision['status']=='PASS_HARMONIZATION_HOLD_REPLICATION_AND_LD'
    cfg=json.loads((ROOT/'config/t17_harmonization.json').read_text())
    assert (cfg['minimum_maf'],cfg['hla_chromosome'],cfg['hla_start'],cfg['hla_end'],
            cfg['bbj_cohort_n_cases'],cfg['bbj_cohort_n_controls'],cfg['exclude_hla_in_both_builds']) == (
                .005,6,25000000,34000000,2176,210277,True), 'Frozen-contract constants changed'
    db=sqlite3.connect(f'file:{out}/harmonization.sqlite?mode=ro',uri=True)
    raw=ROOT/'raw_data/external_gwas/bbj_2019_graves/GD.auto.rsq07.mac10.txt.gz'
    errors=collections.Counter(); counts=collections.Counter()
    cur=iter(db.execute('SELECT id,c37,p37,ref37,alt37,b,se,p,eaf,n,rsq,mac FROM v ORDER BY id'))
    expected=next(cur,None)
    with gzip.open(raw,'rt') as f:
        fields=next(f).split()
        for i,line in enumerate(f,1):
            counts['raw_rows_scanned']+=1
            if expected is None or i!=expected[0]:continue
            r=dict(zip(fields,line.split()))
            observed=(i,int(r['CHR']),int(r['POS']),r['Allele1'],r['Allele2'],float(r['BETA']),float(r['SE']),r['p.value'],float(r['AF_Allele2']),int(r['N']),float(r['Rsq']),float(r['MAC']))
            if observed!=expected:errors['retained_raw_mismatch']+=1
            counts['retained_raw_verified']+=1
            expected=next(cur,None)
    if expected is not None:errors['unmatched_retained_rows']+=1
    print('All retained rows checked against raw BBJ',flush=True)
    query="""SELECT v.id,v.c38,v.p38,v.c37,v.p37,v.ref37,v.alt37,v.b,v.se,v.p,v.eaf,v.n,v.rsq,v.mac,v.strand,v.status,
        f.vid,f.rsid,f.ref,f.alt,f.b,f.se,f.pval,f.eaf
        FROM v JOIN fg f ON f.c=v.c38 AND f.p=v.p38
        WHERE v.status IN ('ALIGNED_DIRECT','ALIGNED_SWAP') ORDER BY v.c38,v.p38"""
    previous=(0,0)
    complement={'A':'T','T':'A','G':'C','C':'G'}
    path=out/'bbj_finngen_graves.aligned.grch38.tsv.gz'
    with gzip.open(path,'rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        for r,x in itertools.zip_longest(reader,db.execute(query)):
            if r is None or x is None:errors['output_row_count']+=1;continue
            (i,c,p,c37,p37,ref,alt,b,se,pval,eaf,n,rsq,mac,strand,status,vid,rsid,fr,fa,fb,fs,fp,fe)=x
            if strand=='-':ref,alt=complement[ref],complement[alt]
            direct=(ref,alt)==(fr,fa);swap=(ref,alt)==(fa,fr)
            if direct==swap:errors['orientation_ambiguous_or_mismatch']+=1
            sign=-1 if swap else 1
            expected_status='ALIGNED_SWAP' if swap else 'ALIGNED_DIRECT'
            if status!=expected_status:errors['sql_orientation_mismatch']+=1
            numeric={'source_row':i,'chromosome':c,'position_grch38':p,'bbj_beta_aligned':sign*b,'bbj_se':se,
                'bbj_eaf_aligned':1-eaf if swap else eaf,'bbj_n_total':n,'bbj_rsq':rsq,'bbj_mac':mac,
                'finngen_beta':fb,'finngen_se':fs,'finngen_eaf':fe,'source_position_grch37':p37,
                'bbj_beta_source':b,'bbj_eaf_source':eaf}
            strings={'discovery_variant_id':vid,'discovery_rsid':rsid,'effect_allele':fa,'other_allele':fr,
                'bbj_p_spa':pval,'finngen_p':fp,'source_ref_grch37':x[5],'source_alt_grch37':x[6],
                'lift_strand':strand,'alignment':status,'bbj_n_cases':'2176' if n==212453 else 'NA',
                'bbj_n_controls':'210277' if n==212453 else 'NA'}
            for k,v in numeric.items():
                if float(r[k])!=v:errors['output_mismatch:'+k]+=1
            for k,v in strings.items():
                if r[k]!=v:errors['output_mismatch:'+k]+=1
            if (c,p)<=previous:errors['output_not_strictly_sorted_unique']+=1
            previous=(c,p)
            counts['aligned_rows_verified']+=1
            counts[expected_status+':'+strand]+=1
            if n!=212453:counts['aligned_reduced_n']+=1
    if counts['raw_rows_scanned']!=8678731:errors['raw_count']+=1
    if counts['retained_raw_verified']!=decision['clean_records']:errors['retained_count']+=1
    if counts['aligned_rows_verified']!=decision['aligned_records']:errors['aligned_count']+=1
    result={'status':'PASS' if not errors else 'FAIL','counts':dict(counts),'errors':dict(errors),
        'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'All retained source rows plus every exported aligned row; no replication tests'}
    (out/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    db.close()
    assert not errors


if __name__=='__main__':main()
