#!/usr/bin/env python3
"""Conservative BBJ Graves cleaning and allele alignment; no statistical tests.

Stages are checkpointed against input hashes. Existing source data and prior task
outputs are immutable. An interrupted stage requires a new --output directory.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import sqlite3
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
COMP = str.maketrans('ACGT', 'TGCA')
SCHEMA = ['id', 'c37', 'p37', 'ref37', 'alt37', 'b', 'se', 'p', 'eaf', 'n', 'rsq', 'mac']


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def table(path, fields, rows):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='', **({'compresslevel': 1} if path.suffix == '.gz' else {})) as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(fields)
        for r in rows:
            w.writerow(['NA' if x is None else x for x in r])


def allele_sign(ref, alt, strand, fg_ref, fg_alt):
    if strand not in ['+', '-']:
        return None
    if len(ref) != 1 or len(alt) != 1 or not set(ref + alt) <= set('ACGT') or ref == alt:
        return None
    if ref + alt in ['AT', 'TA', 'CG', 'GC']:
        return None
    if strand == '-':
        ref, alt = ref.translate(COMP), alt.translate(COMP)
    if (ref, alt) == (fg_ref, fg_alt):
        return 1
    if (ref, alt) == (fg_alt, fg_ref):
        return -1
    return None


def primary_reason(c, pos, ref, alt, eaf, beta, se, p, cfg):
    if c not in range(1, 23) or pos <= 0:
        return 'NON_AUTOSOME_OR_INVALID_POSITION'
    if not all(math.isfinite(v) for v in [eaf, beta, se, p]) or not (0 <= eaf <= 1 and se > 0 and 0 <= p <= 1):
        return 'INVALID_NUMERIC'
    if len(ref) != 1 or len(alt) != 1 or not set(ref + alt) <= set('ACGT') or ref == alt:
        return 'NON_CANONICAL_SNV'
    if c == cfg['hla_chromosome'] and cfg['hla_start'] <= pos <= cfg['hla_end']:
        return 'HLA_SOURCE_BUILD'
    if ref + alt in ['AT', 'TA', 'CG', 'GC']:
        return 'PALINDROMIC_SNV'
    if min(eaf, 1 - eaf) < cfg['minimum_maf']:
        return 'MAF_BELOW_0_005'
    return None


def stream_groups(path, fg=False):
    with gzip.open(path, 'rt') as f:
        h = next(f).split()
        for key, rows in itertools.groupby(
                ((i, dict(zip(h, line.split()))) for i, line in enumerate(f, 1)),
                key=lambda x: (x[1]['chrom' if fg else 'CHR'], x[1]['pos' if fg else 'POS'])):
            yield key, list(rows)


def clean(db, out, cfg):
    db.execute('CREATE TABLE v(id INTEGER PRIMARY KEY,c37 INTEGER,p37 INTEGER,ref37 TEXT,alt37 TEXT,b REAL,se REAL,p TEXT,eaf REAL,n INTEGER,rsq REAL,mac REAL,status TEXT,c38 INTEGER,p38 INTEGER,strand TEXT,ref38 TEXT,alt38 TEXT)')
    counts = collections.Counter(); batch = []
    with gzip.open(out / 'bbj_cleaning_exclusions.tsv.gz', 'wt', compresslevel=1) as audit:
        w = csv.writer(audit, delimiter='\t', lineterminator='\n')
        w.writerow(['source_row', 'source_variant_id', 'chromosome', 'position_grch37', 'reason'])
        for (chrom, position), rows in stream_groups(ROOT / cfg['bbj_input']):
            for i, r in rows:
                counts['raw_rows'] += 1
                c, pos = int(chrom), int(position)
                reason = 'MULTIPLE_RECORDS_AT_SOURCE_POSITION' if len(rows) > 1 else None
                reason = reason or primary_reason(c, pos, r['Allele1'], r['Allele2'],
                    float(r['AF_Allele2']), float(r['BETA']), float(r['SE']), float(r['p.value']), cfg)
                if not reason and (float(r['Rsq']) < cfg['minimum_rsq'] or float(r['MAC']) < cfg['minimum_mac']):
                    reason = 'IMPUTATION_OR_MAC_FILTER'
                if not reason and r['Is.SPA.converge'] != '1':
                    reason = 'SPA_NOT_CONVERGED'
                if not reason and not (0 < int(r['N']) <= 212453):
                    reason = 'INVALID_VARIANT_N'
                if reason:
                    counts[reason] += 1
                    w.writerow([i, r['SNPID'], c, pos, reason])
                else:
                    batch.append((i, c, pos, r['Allele1'], r['Allele2'], float(r['BETA']),
                        float(r['SE']), r['p.value'], float(r['AF_Allele2']), int(r['N']), float(r['Rsq']), float(r['MAC'])))
                    counts['retained'] += 1
                if len(batch) >= 20000:
                    db.executemany('INSERT INTO v(' + ','.join(SCHEMA) + ') VALUES(' + ','.join('?' * len(SCHEMA)) + ')', batch)
                    batch.clear()
                if i % 2000000 == 0:
                    print(f'Clean: {i:,} source rows', flush=True)
    if batch:
        db.executemany('INSERT INTO v(' + ','.join(SCHEMA) + ') VALUES(' + ','.join('?' * len(SCHEMA)) + ')', batch)
    db.execute('CREATE UNIQUE INDEX v_source_position ON v(c37,p37)')
    db.commit()
    table(out / 'bbj_clean.grch37.tsv.gz', SCHEMA,
          db.execute('SELECT ' + ','.join(SCHEMA) + ' FROM v ORDER BY c37,p37'))
    with (out / 'bbj_clean.grch37.bed').open('w') as f:
        for i,c,p in db.execute('SELECT id,c37,p37 FROM v ORDER BY c37,p37'):
            f.write(f'chr{c}\t{p-1}\t{p}\t{i}\t0\t+\n')
    table(out / 'cleaning_counts.tsv', ['reason', 'records'], sorted(counts.items()))
    assert counts['raw_rows'] == counts['retained'] + sum(v for k,v in counts.items() if k not in ['raw_rows', 'retained'])
    assert counts['raw_rows'] == 8678731
    return dict(counts)


def load_bed(db, name, path):
    db.execute(f'CREATE TABLE {name}(id INTEGER,c TEXT,p INTEGER,strand TEXT,width INTEGER)')
    batch = []
    with path.open() as f:
        for line in f:
            c,start,end,i,score,strand = line.rstrip().split('\t')[:6]
            batch.append((int(i),c.removeprefix('chr'),int(start)+1,strand,int(end)-int(start)))
            if len(batch) == 20000:
                db.executemany(f'INSERT INTO {name} VALUES(?,?,?,?,?)', batch);batch.clear()
    if batch:
        db.executemany(f'INSERT INTO {name} VALUES(?,?,?,?,?)', batch)
    db.execute(f'CREATE INDEX {name}_id ON {name}(id)')
    db.execute(f'CREATE TABLE {name}_one AS SELECT id,c,p,strand,width FROM {name} GROUP BY id HAVING COUNT(*)=1')
    db.execute(f'CREATE UNIQUE INDEX {name}_one_id ON {name}_one(id)')
    db.commit()


def lift(db, out, cfg):
    exe = str(ROOT / cfg['liftover_executable'])
    with (out / 'liftover.log').open('w') as log:
        subprocess.run([exe, '-minMatch=1.0', '-multiple', '-noSerial',
            str(out/'bbj_clean.grch37.bed'), str(ROOT/cfg['forward_chain']),
            str(out/'forward.grch38.bed'), str(out/'forward.unmapped.bed')], check=True, stdout=log, stderr=log)
        load_bed(db, 'fw', out/'forward.grch38.bed')
        db.execute("UPDATE v SET status='FORWARD_UNMAPPED'")
        db.execute("UPDATE v SET status='FORWARD_MULTIMAP' WHERE id IN (SELECT id FROM fw GROUP BY id HAVING COUNT(*)>1)")
        db.execute("""UPDATE v SET c38=CAST(f.c AS INT),p38=f.p,strand=f.strand,
             status=CASE WHEN f.c=CAST(v.c37 AS TEXT) AND f.width=1 AND f.strand IN ('+','-')
             THEN 'ROUNDTRIP_PENDING' ELSE 'FORWARD_NONPRIMARY_OR_CHANGED_CHROMOSOME' END
             FROM fw_one f WHERE v.id=f.id""")
        db.commit()
        with (out/'reverse_input.grch38.bed').open('w') as f:
            for i,c,p,s in db.execute("SELECT id,c38,p38,strand FROM v WHERE status='ROUNDTRIP_PENDING'"):
                f.write(f'chr{c}\t{p-1}\t{p}\t{i}\t0\t{s}\n')
        subprocess.run([exe, '-minMatch=1.0', '-multiple', '-noSerial',
            str(out/'reverse_input.grch38.bed'), str(ROOT/cfg['reverse_chain']),
            str(out/'reverse.grch37.bed'), str(out/'reverse.unmapped.bed')], check=True, stdout=log, stderr=log)
    load_bed(db, 'rv', out/'reverse.grch37.bed')
    db.execute("UPDATE v SET status='ROUNDTRIP_UNMAPPED' WHERE status='ROUNDTRIP_PENDING'")
    db.execute("UPDATE v SET status='ROUNDTRIP_MULTIMAP' WHERE id IN (SELECT id FROM rv GROUP BY id HAVING COUNT(*)>1)")
    db.execute("""UPDATE v SET status=CASE WHEN r.c=CAST(v.c37 AS TEXT) AND r.p=v.p37 AND r.width=1 AND r.strand='+'
        THEN 'READY' ELSE 'ROUNDTRIP_MISMATCH' END FROM rv_one r WHERE v.id=r.id""")
    db.execute("UPDATE v SET status='HLA_TARGET_BUILD' WHERE status='READY' AND c38=6 AND p38 BETWEEN ? AND ?", (cfg['hla_start'],cfg['hla_end']))
    db.execute("CREATE TABLE target_collisions AS SELECT c38,p38 FROM v WHERE status='READY' GROUP BY c38,p38 HAVING COUNT(*)>1")
    db.execute("CREATE UNIQUE INDEX target_collisions_pos ON target_collisions(c38,p38)")
    db.execute("UPDATE v SET status='TARGET_POSITION_COLLISION' WHERE status='READY' AND EXISTS(SELECT 1 FROM target_collisions t WHERE t.c38=v.c38 AND t.p38=v.p38)")
    for old,new in [('ref37','ref38'),('alt37','alt38')]:
        db.execute(f"UPDATE v SET {new}=CASE WHEN strand='+' THEN {old} ELSE CASE {old} WHEN 'A' THEN 'T' WHEN 'T' THEN 'A' WHEN 'C' THEN 'G' WHEN 'G' THEN 'C' END END WHERE status='READY'")
    db.execute('CREATE INDEX v_target_pos ON v(c38,p38)')
    db.commit()
    counts = dict(db.execute('SELECT status,COUNT(*) FROM v GROUP BY status'))
    table(out/'lift_counts.tsv', ['status','records'], sorted(counts.items()))
    return counts


def align(db, out, cfg):
    db.execute('CREATE TABLE fg(c INTEGER,p INTEGER,ref TEXT,alt TEXT,b REAL,se REAL,pval TEXT,eaf REAL,vid TEXT,rsid TEXT,PRIMARY KEY(c,p))')
    db.execute('CREATE TABLE fg_reject(c INTEGER,p INTEGER,reason TEXT,PRIMARY KEY(c,p))')
    counts=collections.Counter();batch=[];reject=[]
    for (chrom,position),rows in stream_groups(ROOT/cfg['discovery_input'], fg=True):
        c,pos=int(chrom),int(position)
        counts['discovery_source_rows']+=len(rows)
        i,r=rows[0]
        reason='MULTIPLE_RECORDS_AT_DISCOVERY_POSITION' if len(rows)>1 else primary_reason(
            c,pos,r['ref'],r['alt'],float(r['eaf']),float(r['beta']),float(r['se']),float(r['p']),cfg)
        if reason:
            reject.append((c,pos,reason));counts['discovery_excluded_rows']+=len(rows)
        else:
            batch.append((c,pos,r['ref'],r['alt'],float(r['beta']),float(r['se']),r['p'],float(r['eaf']),r['variant_id'],r['rsid']))
            counts['discovery_eligible_rows']+=1
        if len(batch)+len(reject)>=20000:
            db.executemany('INSERT INTO fg VALUES(?,?,?,?,?,?,?,?,?,?)',batch);batch.clear()
            db.executemany('INSERT INTO fg_reject VALUES(?,?,?)',reject);reject.clear()
        if i%2000000==0:print(f'Align: {i:,} discovery rows',flush=True)
    if batch:db.executemany('INSERT INTO fg VALUES(?,?,?,?,?,?,?,?,?,?)',batch)
    if reject:db.executemany('INSERT INTO fg_reject VALUES(?,?,?)',reject)
    db.execute('CREATE UNIQUE INDEX fg_variant_id ON fg(vid)')
    assert db.execute('SELECT COUNT(*) FROM fg JOIN fg_reject USING(c,p)').fetchone()[0]==0
    db.execute("""UPDATE v SET status=CASE WHEN EXISTS(SELECT 1 FROM fg_reject f WHERE f.c=v.c38 AND f.p=v.p38)
        THEN 'DISCOVERY_POSITION_EXCLUDED' ELSE 'NO_DISCOVERY_POSITION' END
        WHERE status='READY' AND NOT EXISTS(SELECT 1 FROM fg f WHERE f.c=v.c38 AND f.p=v.p38)""")
    db.execute("""UPDATE v SET status=CASE WHEN v.ref38=f.ref AND v.alt38=f.alt THEN 'ALIGNED_DIRECT'
        WHEN v.ref38=f.alt AND v.alt38=f.ref THEN 'ALIGNED_SWAP' ELSE 'ALLELE_MISMATCH' END
        FROM fg f WHERE status='READY' AND f.c=v.c38 AND f.p=v.p38""")
    db.commit()
    cols=['source_row','chromosome','position_grch38','discovery_variant_id','discovery_rsid',
          'effect_allele','other_allele','bbj_beta_aligned','bbj_se','bbj_p_spa','bbj_eaf_aligned',
          'bbj_n_total','bbj_n_cases','bbj_n_controls','bbj_rsq','bbj_mac','finngen_beta','finngen_se',
          'finngen_p','finngen_eaf','source_position_grch37','source_ref_grch37','source_alt_grch37',
          'lift_strand','alignment','bbj_beta_source','bbj_eaf_source']
    query="""SELECT v.id,v.c38,v.p38,f.vid,f.rsid,f.alt,f.ref,
        CASE WHEN v.status='ALIGNED_SWAP' THEN -v.b ELSE v.b END,v.se,v.p,
        CASE WHEN v.status='ALIGNED_SWAP' THEN 1-v.eaf ELSE v.eaf END,v.n,
        CASE WHEN v.n=212453 THEN 2176 END,CASE WHEN v.n=212453 THEN 210277 END,
        v.rsq,v.mac,f.b,f.se,f.pval,f.eaf,v.p37,v.ref37,v.alt37,v.strand,v.status,v.b,v.eaf
        FROM v JOIN fg f ON v.c38=f.c AND v.p38=f.p WHERE v.status IN ('ALIGNED_DIRECT','ALIGNED_SWAP')
        ORDER BY v.c38,v.p38"""
    table(out/'bbj_finngen_graves.aligned.grch38.tsv.gz',cols,db.execute(query))
    table(out/'bbj_mapping_alignment_audit.tsv.gz',
          ['source_row','chromosome_grch37','position_grch37','ref_grch37','alt_grch37','chromosome_grch38','position_grch38','lift_strand','status'],
          db.execute('SELECT id,c37,p37,ref37,alt37,c38,p38,strand,status FROM v ORDER BY c37,p37'))
    counts.update(dict(db.execute('SELECT status,COUNT(*) FROM v GROUP BY status')))
    table(out/'alignment_counts.tsv',['status','records'],sorted(counts.items()))
    return dict(counts)


def finish(db,out,cfg,provenance,state):
    rows=[]
    with (ROOT/'outcome_data/t07/genetic_signal_five_classification.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            match=db.execute("SELECT v.status,v.id FROM v JOIN fg f ON f.c=v.c38 AND f.p=v.p38 WHERE f.vid=? AND v.status IN ('ALIGNED_DIRECT','ALIGNED_SWAP')",(r['representative_variant'],)).fetchone()
            rows.append([r['signal_unit_id'],r['five_class_primary'],r['representative_variant'],
                         'IN_ALIGNED_PRIMARY_SET' if match else 'NOT_IN_ALIGNED_PRIMARY_SET',match[1] if match else None,
                         'Availability only; absence may reflect filters/mapping/ancestry, not non-replication'])
    table(out/'t07_representative_variant_coverage.tsv',
          ['signal_unit_id','five_class_primary','representative_variant','coverage','bbj_source_row','interpretation'],rows)
    checks=[]
    def check(name,ok):checks.append({'check':name,'passed':bool(ok)})
    check('raw_partition',state['clean']['raw_rows']==8678731)
    check('clean_unique_source_positions',db.execute('SELECT COUNT(*) FROM v').fetchone()[0]==state['clean']['retained'])
    check('mapping_partition',sum(state['lift'].values())==state['clean']['retained'])
    check('final_partition',sum(n for s,n in db.execute('SELECT status,COUNT(*) FROM v GROUP BY status'))==state['clean']['retained'])
    check('no_pending_states',db.execute("SELECT COUNT(*) FROM v WHERE status IN ('READY','ROUNDTRIP_PENDING') OR status IS NULL").fetchone()[0]==0)
    check('aligned_nonempty',db.execute("SELECT COUNT(*) FROM v WHERE status LIKE 'ALIGNED_%'").fetchone()[0]>0)
    check('aligned_unique_target_positions',db.execute("SELECT COUNT(*) FROM (SELECT c38,p38 FROM v WHERE status LIKE 'ALIGNED_%' GROUP BY c38,p38 HAVING COUNT(*)>1)").fetchone()[0]==0)
    check('aligned_exact_roundtrip',db.execute("SELECT COUNT(*) FROM v LEFT JOIN rv_one r ON v.id=r.id WHERE v.status LIKE 'ALIGNED_%' AND (r.id IS NULL OR r.c!=CAST(v.c37 AS TEXT) OR r.p!=v.p37 OR r.strand!='+' OR r.width!=1)").fetchone()[0]==0)
    check('no_hla_in_either_build',db.execute("SELECT COUNT(*) FROM v WHERE status LIKE 'ALIGNED_%' AND ((c37=6 AND p37 BETWEEN 25000000 AND 34000000) OR (c38=6 AND p38 BETWEEN 25000000 AND 34000000))").fetchone()[0]==0)
    check('no_palindromic_in_aligned',db.execute("SELECT COUNT(*) FROM v WHERE status LIKE 'ALIGNED_%' AND ref37||alt37 IN ('AT','TA','CG','GC')").fetchone()[0]==0)
    check('t07_coverage_complete',len(rows)==304)
    for p,expected in provenance.items():check('input_unchanged:'+p,digest(ROOT/p)==expected)
    (out/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    statuses=dict(db.execute('SELECT status,COUNT(*) FROM v GROUP BY status'))
    decision={'completed_at_utc':datetime.now(timezone.utc).isoformat(),
        'status':'PASS_HARMONIZATION_HOLD_REPLICATION_AND_LD' if all(x['passed'] for x in checks) else 'HOLD_VALIDATION',
        'source_records':state['clean']['raw_rows'],'clean_records':state['clean']['retained'],
        'aligned_records':sum(v for k,v in statuses.items() if k.startswith('ALIGNED_')),
        'alignment_statuses':statuses,'checks_passed':sum(x['passed'] for x in checks),'checks_total':len(checks),
        'liftover_executed':True,'roundtrip_verified':True,'replication_executed':False,'ld_used':False,
        'holds':['No EAS/Japanese LD or signal-level colocalisation','No replication/heterogeneity/power tests',
                 'Indels, symbolic variants, palindromic SNPs and HLA outside primary set',
                 'Exact phenotype equivalence not established','Independent EUR Graves not acquired',
                 'No independent reference-FASTA validation; accepted pairs require discovery allele concordance']}
    (out/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
    print(json.dumps(decision,indent=2),flush=True)
    assert all(x['passed'] for x in checks)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='outcome_data/t17/harmonization')
    args=ap.parse_args();out=(ROOT/args.output).resolve();out.mkdir(parents=True,exist_ok=True)
    assert out.is_relative_to(ROOT/'outcome_data/t17')
    cfgpath=ROOT/'config/t17_harmonization.json';cfg=json.loads(cfgpath.read_text())
    inputs=[cfgpath,Path(__file__).resolve(),ROOT/'config/analysis.yml',ROOT/'config/traits.tsv',
            ROOT/'outcome_data/t07/genetic_signal_five_classification.tsv',
            *[ROOT/cfg[k] for k in ['bbj_input','discovery_input','forward_chain','reverse_chain','liftover_executable']]]
    provenance={str(p.relative_to(ROOT)):digest(p) for p in inputs}
    statepath=out/'checkpoint.json'
    state=json.loads(statepath.read_text()) if statepath.exists() else {'provenance':provenance}
    assert state['provenance']==provenance,'Inputs/code changed; use a new --output directory'
    if state.get('running'):raise RuntimeError('Interrupted stage: use a new --output directory')
    db=sqlite3.connect(out/'harmonization.sqlite');db.execute('PRAGMA cache_size=-262144')
    db.execute('PRAGMA journal_mode=WAL');db.execute('PRAGMA synchronous=NORMAL')
    for name,fn in [('clean',clean),('lift',lift),('align',align)]:
        if name in state:continue
        state['running']=name;statepath.write_text(json.dumps(state,indent=2)+'\n')
        state[name]=fn(db,out,cfg);state.pop('running');statepath.write_text(json.dumps(state,indent=2)+'\n')
    finish(db,out,cfg,provenance,state)
    table(out/'provenance.tsv',['path','sha256'],provenance.items())
    db.execute('PRAGMA wal_checkpoint(TRUNCATE)');db.close()


if __name__=='__main__':main()
