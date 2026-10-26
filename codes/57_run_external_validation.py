#!/usr/bin/env python3
"""Frozen priority-variant cross-ancestry comparison; never signal replication."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
NORMAL=NormalDist()
STAT_FIELDS=['direction','source_p_bonferroni','normal_p','external_or','external_or_lower95','external_or_upper95',
 'discovery_or','discovery_or_lower95','discovery_or_upper95','difference_beta','difference_se','difference_ci_lower95',
 'difference_ci_upper95','heterogeneity_q','heterogeneity_p','heterogeneity_p_bonferroni','i2_descriptive',
 'equivalence_p','equivalence_p_bonferroni','equivalence_flag','equivalence_beta_lower995','equivalence_beta_upper995',
 'precision_power_or1_10','precision_power_half_discovery','precision_power_full_discovery','minimum_detectable_or80',
 'classification']


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def sf(z):return .5*math.erfc(z/math.sqrt(2))


def precision_power(effect,se,alpha):
    z=NORMAL.inv_cdf(1-alpha/2); delta=abs(effect)/se
    return sf(z-delta)+sf(z+delta)


def mde(se,alpha):
    lo,hi=0.,se*10
    for _ in range(80):
        mid=(lo+hi)/2
        if precision_power(mid,se,alpha)<.8:lo=mid
        else:hi=mid
    return math.exp((lo+hi)/2)


def calculate(b,se,p,bf,sf_,cfg):
    assert all(math.isfinite(x) for x in [b,se,p,bf,sf_]) and se>0 and sf_>0 and 0<=p<=1
    n=cfg['planned_tests']; alpha=cfg['family_alpha']/n; margin=math.log(cfg['equivalence_or_upper'])
    z95=NORMAL.inv_cdf(.975); zeq=NORMAL.inv_cdf(1-alpha)
    difference=b-bf; dse=math.hypot(se,sf_); q=(difference/dse)**2
    hp=math.erfc(math.sqrt(q/2)); ep=max(sf((b+margin)/se),sf((margin-b)/se))
    direction='same' if b*bf>0 else 'opposite' if b*bf<0 else 'zero'
    pp=precision_power(margin,se,alpha)
    equivalent=ep<=alpha
    if p<=alpha and direction=='same':category='same_direction_association_support'
    elif p<=alpha and direction=='opposite':category='opposite_direction_association_support'
    elif equivalent:category='no_primary_support_small_effect_equivalence'
    elif pp<cfg['precision_power_threshold']:category='inconclusive_low_precision_OR1_10_scenario'
    else:category='inconclusive'
    values=[direction,min(1,p*n),math.erfc(abs(b/se)/math.sqrt(2)),math.exp(b),math.exp(b-z95*se),math.exp(b+z95*se),
        math.exp(bf),math.exp(bf-z95*sf_),math.exp(bf+z95*sf_),difference,dse,difference-z95*dse,difference+z95*dse,
        q,hp,min(1,hp*n),max(0,(q-1)/q) if q else 0,ep,min(1,ep*n),equivalent,b-zeq*se,b+zeq*se,pp,
        precision_power(.5*bf,se,alpha),precision_power(bf,se,alpha),mde(se,alpha),category]
    return dict(zip(STAT_FIELDS,values))


def read_targets(path,wanted,key):
    found={}
    with gzip.open(path,'rt') as f:
        fields=next(f).rstrip('\n').split('\t'); index=fields.index(key)
        for line in f:
            parts=line.rstrip('\n').split('\t', index+1)
            if parts[index] not in wanted:continue
            r=dict(zip(fields,line.rstrip('\n').split('\t')))
            assert r[key] not in found,'Duplicate target in '+str(path)
            found[r[key]]=r
    return found


def write_table(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:w.writerow({k:'NA' if r.get(k) is None else r[k] for k in fields})


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',default='outcome_data/external_validation/plan_20260904')
    ap.add_argument('--output',default='outcome_data/external_validation/run_20260904');args=ap.parse_args()
    plan=ROOT/args.plan;out=ROOT/args.output
    assert out.resolve().is_relative_to(ROOT/'outcome_data/external_validation') and not out.exists()
    cfg=json.loads((plan/'plan.json').read_text());frozen=json.loads((plan/'freeze.json').read_text())
    for name,expected in frozen['sha256'].items():assert sha(ROOT/name)==expected,'Frozen input changed: '+name
    assert cfg['method_review']=='PASS_WITH_EXPLICIT_APPROXIMATION_LIMITS'
    with (plan/'targets.tsv').open() as f:targets=list(csv.DictReader(f,delimiter='\t'))
    assert len(targets)==10 and cfg['planned_tests']==20
    wanted={r['representative_variant'] for r in targets};assert len(wanted)==10
    rows=[];provenance=[];checks=[]
    for trait,s in cfg['traits'].items():
        for key in ['decision','independent_audit']:
            d=json.loads((ROOT/s[key]).read_text())
            if key=='decision':assert d['status'].startswith('PASS'),s[key]
            else:assert d['status']=='PASS',s[key]
        print('Extracting frozen targets: '+trait,flush=True)
        fg=read_targets(ROOT/s['discovery_input'],wanted,'variant_id')
        ext=read_targets(ROOT/s['aligned_input'],wanted,'discovery_variant_id')
        assert set(fg)==wanted,'Missing discovery targets'
        for t in targets:
            vid=t['representative_variant'];f=fg[vid];e=ext.get(vid)
            assert (f['ref'],f['alt'])==(t['ref'],t['alt']) and f['effect_allele']==f['alt']
            r={k:None for k in STAT_FIELDS}
            r.update(signal_unit_id=t['signal_unit_id'],five_class_discovery=t['five_class_primary'],trait=trait,
                variant_id=vid,nearest_genes=t['representative_nearest_genes'],effect_allele=f['alt'],other_allele=f['ref'],
                discovery_beta=float(f['beta']),discovery_se=float(f['se']),discovery_p=f['p'],discovery_eaf=float(f['eaf']),
                external_beta=None,external_se=None,external_source_p=None,external_eaf=None,external_n_total=None,
                external_source_row=None,external_alignment=None,external_lift_strand=None,
                availability='not_in_aligned_primary_set',classification='untestable_missing_aligned_target',
                discovery_class_external_status='NOT_REVALIDATED_NO_EXTERNAL_SIGNAL_LEVEL_EVIDENCE')
            if e is not None:
                assert (e['effect_allele'],e['other_allele'])==(f['alt'],f['ref'])
                assert math.isclose(float(e['finngen_beta']),float(f['beta']),rel_tol=1e-12,abs_tol=1e-12)
                assert math.isclose(float(e['finngen_se']),float(f['se']),rel_tol=1e-12,abs_tol=1e-12)
                b,se,p=map(float,[e['bbj_beta_aligned'],e['bbj_se'],e[s['p_column']]])
                r.update(calculate(b,se,p,float(f['beta']),float(f['se']),cfg))
                r.update(external_beta=b,external_se=se,external_source_p=e[s['p_column']],external_eaf=float(e['bbj_eaf_aligned']),
                    external_n_total=e.get('bbj_n_total','NA'),external_source_row=e['source_row'],external_alignment=e['alignment'],
                    external_lift_strand=e['lift_strand'],availability='aligned_exact_target')
            rows.append(r)
        for key in ['discovery_input','aligned_input','decision','independent_audit']:
            path=ROOT/s[key];provenance.append({'path':s[key],'sha256':sha(path)})
    assert len(rows)==20 and len({(r['signal_unit_id'],r['trait']) for r in rows})==20
    counts=collections.Counter((r['trait'],r['classification']) for r in rows)
    for name,expected in frozen['sha256'].items():assert sha(ROOT/name)==expected,'Frozen input changed during run'
    out.mkdir(parents=True)
    write_table(out/'priority_variant_validation.tsv',rows,list(rows[0]))
    write_table(out/'classification_counts.tsv',[{'trait':t,'classification':c,'n':n} for (t,c),n in sorted(counts.items())],['trait','classification','n'])
    provenance.append({'path':str(Path(__file__).resolve().relative_to(ROOT)),'sha256':sha(Path(__file__))})
    write_table(out/'provenance.tsv',provenance,['path','sha256'])
    decision={'status':'PASS_FEASIBLE_VARIANT_VALIDATION_HOLD_SIGNAL_LEVEL_AND_T18_SYNTHESIS',
        'completed_at_utc':datetime.now(timezone.utc).isoformat(),'planned_tests':20,
        'testable_tests':sum(r['availability']=='aligned_exact_target' for r in rows),
        'counts':[{'trait':t,'classification':c,'n':n} for (t,c),n in sorted(counts.items())],
        'source_p_association_executed':True,'normal_approximation_heterogeneity_executed':True,
        'normal_approximation_equivalence_executed':True,'signal_level_replication_executed':False,
        'cross_disease_difference_executed':False,'pooled_meta_analysis_executed':False,
        'limitations':['Normal CI/Q/TOST/SE precision scenarios do not inherit SPA calibration',
            'No validated EAS/Japanese LD; no signal-level colocalisation',
            'Phenotype equivalence not established; cross-ancestry effects can differ via LD or biology',
            'BBJ disease cohorts overlap; no cross-disease contrast without covariance',
            'OR1.10 precision and equivalence are prespecified operational margins, not biological absence',
            'Discovery effect scenarios may be optimistic from winner selection; no external observed-power claim',
            'Independent EUR sources unavailable; original T07 classes are not revalidated']}
    (out/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
    print(json.dumps(decision,indent=2),flush=True)


if __name__=='__main__':main()
