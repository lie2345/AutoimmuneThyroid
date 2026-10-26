#!/usr/bin/env python3
"""Conservative BBJ Hashimoto (T16) GRCh37-to-GRCh38 harmonisation.

The BBJ EA/OA fields are study effect/other alleles.  They are deliberately not
renamed REF/ALT: no reference-genome claim is available from this input.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import math
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMP = str.maketrans('ACGT', 'TGCA')


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_table(path, fields, rows):
    opener = gzip.open if path.suffix == '.gz' else open
    kwargs = {'newline': ''}
    if path.suffix == '.gz':
        kwargs['compresslevel'] = 1
    with opener(path, 'wt', **kwargs) as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        for row in rows:
            writer.writerow(['NA' if value is None else value for value in row])


def parse_chromosome(value):
    value = value.removeprefix('chr')
    return int(value) if value.isdigit() else None


def is_snv(a, b):
    return len(a) == len(b) == 1 and a != b and set(a + b) <= set('ACGT')


def is_palindromic(a, b):
    return a + b in {'AT', 'TA', 'CG', 'GC'}


def primary_reason(c, pos, ea, oa, eaf, beta, se, pval, cfg, build):
    if c not in range(1, 23) or pos <= 0:
        return 'NON_AUTOSOME_OR_INVALID_POSITION'
    if not all(math.isfinite(x) for x in (eaf, beta, se, pval)) or not (0 <= eaf <= 1 and se > 0 and 0 <= pval <= 1):
        return 'INVALID_NUMERIC'
    if not is_snv(ea, oa):
        return 'NON_CANONICAL_SNV'
    if c == cfg['hla_chromosome'] and cfg['hla_start'] <= pos <= cfg['hla_end']:
        return 'HLA_' + build.upper()
    if is_palindromic(ea, oa):
        return 'PALINDROMIC_SNV'
    if min(eaf, 1 - eaf) < cfg['minimum_maf']:
        return 'MAF_BELOW_0_005'
    return None


def source_order(path):
    """Return row count, X count, and whether `(chromosome, position)` is numeric-sorted."""
    rows = x_rows = 0
    ordered, previous = True, None
    with gzip.open(path, 'rt') as handle:
        fields = next(handle).rstrip('\n').split('\t')
        for line in handle:
            rows += 1
            record = dict(zip(fields, line.rstrip('\n').split('\t')))
            raw_c = record['chromosome']
            if raw_c.upper().removeprefix('CHR') == 'X':
                x_rows += 1
            c = parse_chromosome(raw_c)
            try:
                pos = int(record['base_pair_location'])
            except ValueError:
                ordered = False
                continue
            key = (c if c is not None else 99, pos)
            if previous is not None and key < previous:
                ordered = False
            previous = key
    return rows, x_rows, ordered


def load_bed(db, name, path):
    db.execute(f'CREATE TABLE {name}(id INTEGER,c TEXT,p INTEGER,strand TEXT,width INTEGER)')
    batch = []
    with path.open() as handle:
        for line in handle:
            c, start, end, ident, _, strand = line.rstrip('\n').split('\t')[:6]
            batch.append((int(ident), c.removeprefix('chr'), int(start) + 1, strand, int(end) - int(start)))
            if len(batch) == 20000:
                db.executemany(f'INSERT INTO {name} VALUES(?,?,?,?,?)', batch); batch.clear()
    if batch:
        db.executemany(f'INSERT INTO {name} VALUES(?,?,?,?,?)', batch)
    db.execute(f'CREATE INDEX {name}_id ON {name}(id)')
    db.execute(f'CREATE TABLE {name}_one AS SELECT id,c,p,strand,width FROM {name} GROUP BY id HAVING COUNT(*)=1')
    db.execute(f'CREATE UNIQUE INDEX {name}_one_id ON {name}_one(id)')
    db.commit()


def prepare_source(db, out, cfg):
    raw = ROOT / cfg['bbj_input']
    rows, x_rows, ordered = source_order(raw)
    if rows != cfg['bbj_expected_rows'] or x_rows != cfg['bbj_expected_x_rows']:
        raise RuntimeError(f'Unexpected source row counts: {rows=}, {x_rows=}')
    db.execute('''CREATE TABLE v(id INTEGER PRIMARY KEY,c37 INTEGER,p37 INTEGER,ea37 TEXT,oa37 TEXT,
                  beta REAL,se REAL,p TEXT,eaf REAL,status TEXT,c38 INTEGER,p38 INTEGER,strand TEXT,ea38 TEXT,oa38 TEXT)''')
    counts, batch = collections.Counter(), []
    with gzip.open(out / 'bbj_cleaning_exclusions.tsv.gz', 'wt', compresslevel=1, newline='') as audit, gzip.open(raw, 'rt') as handle:
        writer = csv.writer(audit, delimiter='\t', lineterminator='\n')
        writer.writerow(['source_row', 'source_variant_id', 'chromosome', 'position_grch37', 'reason'])
        fields = next(handle).rstrip('\n').split('\t')
        # Numeric order was proved in a complete first pass.  Adjacent groups are
        # therefore full position groups, rather than an unverified groupby guess.
        if not ordered:
            raise RuntimeError('Source position order is not numeric-sorted; this run refuses unsafe adjacent grouping')
        parsed = ((source_row, dict(zip(fields, line.rstrip('\n').split('\t')))) for source_row, line in enumerate(handle, 1))
        def pos_key(item):
            r = item[1]
            return (parse_chromosome(r['chromosome']), int(r['base_pair_location']))
        for (c, pos), group in itertools.groupby(parsed, key=pos_key):
            group = list(group)
            multiple = len(group) > 1
            for source_row, record in group:
                counts['raw_rows'] += 1
                try:
                    reason = 'MULTIPLE_RECORDS_AT_SOURCE_POSITION' if multiple else primary_reason(c, pos, record['effect_allele'], record['other_allele'],
                        float(record['effect_allele_frequency']), float(record['beta']), float(record['standard_error']), float(record['p_value']), cfg, 'source')
                except (ValueError, TypeError):
                    reason = 'INVALID_NUMERIC'
                if reason:
                    counts[reason] += 1
                    writer.writerow([source_row, record['variant_id'], record['chromosome'], record['base_pair_location'], reason])
                else:
                    batch.append((source_row, c, pos, record['effect_allele'], record['other_allele'], float(record['beta']), float(record['standard_error']), record['p_value'], float(record['effect_allele_frequency']), 'CANDIDATE'))
                if len(batch) >= 20000:
                    db.executemany('INSERT INTO v(id,c37,p37,ea37,oa37,beta,se,p,eaf,status) VALUES(?,?,?,?,?,?,?,?,?,?)', batch); batch.clear()
            if counts['raw_rows'] % 2000000 < len(group): print(f'Source pass: {counts["raw_rows"]:,}', flush=True)
    if batch: db.executemany('INSERT INTO v(id,c37,p37,ea37,oa37,beta,se,p,eaf,status) VALUES(?,?,?,?,?,?,?,?,?,?)', batch)
    counts['candidate_before_source_duplicates'] = counts['retained'] = db.execute("SELECT COUNT(*) FROM v WHERE status='CANDIDATE'").fetchone()[0]
    db.execute('CREATE UNIQUE INDEX v_source_position ON v(c37,p37)')
    db.commit()
    write_table(out / 'source_order_check.tsv', ['input', 'rows', 'x_rows', 'numeric_chromosome_position_sorted'],
                [(cfg['bbj_input'], rows, x_rows, str(ordered).lower())])
    write_table(out / 'bbj_clean.grch37.tsv.gz', ['source_row','chromosome','position_grch37','source_effect_allele','source_other_allele','bbj_beta_source','bbj_se','bbj_p','bbj_eaf_source'],
                db.execute("SELECT id,c37,p37,ea37,oa37,beta,se,p,eaf FROM v WHERE status='CANDIDATE' ORDER BY c37,p37"))
    with (out / 'bbj_clean.grch37.bed').open('w') as handle:
        for ident, c, pos in db.execute("SELECT id,c37,p37 FROM v WHERE status='CANDIDATE' ORDER BY c37,p37"):
            handle.write(f'chr{c}\t{pos-1}\t{pos}\t{ident}\t0\t+\n')
    write_table(out / 'cleaning_counts.tsv', ['reason','records'], sorted(counts.items()))
    counts['source_x_rows'] = x_rows
    return dict(counts)


def lift(db, out, cfg):
    exe = str(ROOT / cfg['liftover_executable'])
    with (out / 'liftover.log').open('w') as log:
        subprocess.run([exe, '-minMatch=1.0', '-multiple', '-noSerial', str(out/'bbj_clean.grch37.bed'),
            str(ROOT/cfg['forward_chain']), str(out/'forward.grch38.bed'), str(out/'forward.unmapped.bed')], check=True, stdout=log, stderr=log)
        load_bed(db, 'fw', out/'forward.grch38.bed')
        db.execute("UPDATE v SET status='FORWARD_UNMAPPED' WHERE status='CANDIDATE'")
        db.execute("UPDATE v SET status='FORWARD_MULTIMAP' WHERE status='FORWARD_UNMAPPED' AND id IN (SELECT id FROM fw GROUP BY id HAVING COUNT(*)>1)")
        db.execute("""UPDATE v SET c38=CAST(f.c AS INTEGER),p38=f.p,strand=f.strand,
          status=CASE WHEN f.c=CAST(v.c37 AS TEXT) AND f.width=1 AND f.strand IN ('+','-') THEN 'ROUNDTRIP_PENDING'
          ELSE 'FORWARD_NONPRIMARY_OR_CHANGED_CHROMOSOME' END FROM fw_one f WHERE v.status='FORWARD_UNMAPPED' AND v.id=f.id""")
        db.commit()
        with (out / 'reverse_input.grch38.bed').open('w') as handle:
            for ident, c, pos, strand in db.execute("SELECT id,c38,p38,strand FROM v WHERE status='ROUNDTRIP_PENDING'"):
                handle.write(f'chr{c}\t{pos-1}\t{pos}\t{ident}\t0\t{strand}\n')
        subprocess.run([exe, '-minMatch=1.0', '-multiple', '-noSerial', str(out/'reverse_input.grch38.bed'),
            str(ROOT/cfg['reverse_chain']), str(out/'reverse.grch37.bed'), str(out/'reverse.unmapped.bed')], check=True, stdout=log, stderr=log)
    load_bed(db, 'rv', out/'reverse.grch37.bed')
    db.execute("UPDATE v SET status='ROUNDTRIP_UNMAPPED' WHERE status='ROUNDTRIP_PENDING'")
    db.execute("UPDATE v SET status='ROUNDTRIP_MULTIMAP' WHERE status='ROUNDTRIP_UNMAPPED' AND id IN (SELECT id FROM rv GROUP BY id HAVING COUNT(*)>1)")
    db.execute("""UPDATE v SET status=CASE WHEN r.c=CAST(v.c37 AS TEXT) AND r.p=v.p37 AND r.width=1 AND r.strand='+' THEN 'READY'
        ELSE 'ROUNDTRIP_MISMATCH' END FROM rv_one r WHERE v.status='ROUNDTRIP_UNMAPPED' AND v.id=r.id""")
    db.execute("UPDATE v SET status='HLA_TARGET_BUILD' WHERE status='READY' AND c38=6 AND p38 BETWEEN ? AND ?", (cfg['hla_start'],cfg['hla_end']))
    db.execute("CREATE TABLE target_collisions AS SELECT c38,p38 FROM v WHERE status='READY' GROUP BY c38,p38 HAVING COUNT(*)>1")
    db.execute("CREATE UNIQUE INDEX target_collisions_position ON target_collisions(c38,p38)")
    db.execute("UPDATE v SET status='TARGET_POSITION_COLLISION' WHERE status='READY' AND EXISTS(SELECT 1 FROM target_collisions t WHERE t.c38=v.c38 AND t.p38=v.p38)")
    db.execute("""UPDATE v SET ea38=CASE WHEN strand='+' THEN ea37 ELSE CASE ea37 WHEN 'A' THEN 'T' WHEN 'T' THEN 'A' WHEN 'C' THEN 'G' WHEN 'G' THEN 'C' END END,
        oa38=CASE WHEN strand='+' THEN oa37 ELSE CASE oa37 WHEN 'A' THEN 'T' WHEN 'T' THEN 'A' WHEN 'C' THEN 'G' WHEN 'G' THEN 'C' END END WHERE status='READY'""")
    db.execute("CREATE INDEX v_target_position ON v(c38,p38)")
    db.commit()
    counts = dict(db.execute('SELECT status,COUNT(*) FROM v GROUP BY status'))
    write_table(out / 'lift_counts.tsv', ['status','records'], sorted(counts.items()))
    return counts


def load_discovery(db, out, cfg):
    """Load target rows first, then determine duplicates globally (input is not assumed sorted)."""
    db.execute('''CREATE TABLE fg_raw(row_id INTEGER PRIMARY KEY,c INTEGER,p INTEGER,ref TEXT,alt TEXT,beta REAL,se REAL,pval TEXT,eaf REAL,vid TEXT)''')
    db.execute('CREATE TABLE fg_pos(c INTEGER,p INTEGER,row_count INTEGER,PRIMARY KEY(c,p))')
    counts, batch, positions = collections.Counter(), [], []
    with gzip.open(ROOT/cfg['discovery_input'], 'rt') as handle:
        fields = next(handle).rstrip('\n').split('\t')
        for row_id, line in enumerate(handle, 1):
            record = dict(zip(fields, line.rstrip('\n').split('\t')))
            counts['discovery_source_rows'] += 1
            c, pos = parse_chromosome(record['chrom']), None
            try: pos = int(record['pos'])
            except ValueError: pass
            if c is not None and pos is not None: positions.append((c, pos))
            try:
                batch.append((row_id, c, pos, record['ref'], record['alt'], float(record['beta']), float(record['se']), record['p'], float(record['eaf']), record['variant_id']))
            except (ValueError, TypeError):
                # Retain a syntactically valid row with NaNs so primary_reason records INVALID_NUMERIC below.
                batch.append((row_id, c, pos, record['ref'], record['alt'], float('nan'), float('nan'), record['p'], float('nan'), record['variant_id']))
            if len(batch) >= 20000:
                db.executemany('INSERT INTO fg_raw VALUES(?,?,?,?,?,?,?,?,?,?)', batch); batch.clear()
            if len(positions) >= 20000:
                db.executemany('''INSERT INTO fg_pos(c,p,row_count) VALUES(?,?,1) ON CONFLICT(c,p) DO UPDATE SET row_count=row_count+1''', positions); positions.clear()
            if row_id % 2000000 == 0: print(f'Discovery load: {row_id:,}', flush=True)
    if batch: db.executemany('INSERT INTO fg_raw VALUES(?,?,?,?,?,?,?,?,?,?)', batch)
    if positions: db.executemany('''INSERT INTO fg_pos(c,p,row_count) VALUES(?,?,1) ON CONFLICT(c,p) DO UPDATE SET row_count=row_count+1''', positions)
    db.execute('''CREATE TABLE fg AS SELECT r.c,r.p,r.ref,r.alt,r.beta,r.se,r.pval,r.eaf,r.vid
        FROM fg_raw r JOIN fg_pos q ON q.c=r.c AND q.p=r.p WHERE q.row_count=1''')
    db.execute('CREATE UNIQUE INDEX fg_position ON fg(c,p)')
    db.execute('CREATE TABLE fg_reject(c INTEGER,p INTEGER,reason TEXT,PRIMARY KEY(c,p))')
    insert = []
    for c,p,ref,alt,beta,se,pval,eaf in db.execute('SELECT c,p,ref,alt,beta,se,pval,eaf FROM fg'):
        try: reason = primary_reason(c,p,ref,alt,eaf,beta,se,float(pval),cfg,'target')
        except (ValueError, TypeError): reason = 'INVALID_NUMERIC'
        if reason: insert.append((c,p,reason))
        if len(insert) >= 20000:
            db.executemany('INSERT INTO fg_reject VALUES(?,?,?)', insert); insert.clear()
    if insert: db.executemany('INSERT INTO fg_reject VALUES(?,?,?)', insert)
    duplicate_target_rows = db.execute('SELECT COALESCE(SUM(row_count),0) FROM fg_pos WHERE row_count>1').fetchone()[0]
    counts['MULTIPLE_RECORDS_AT_DISCOVERY_POSITION'] = duplicate_target_rows
    counts['discovery_unique_positions'] = db.execute('SELECT COUNT(*) FROM fg').fetchone()[0]
    counts['discovery_eligible_positions'] = db.execute('SELECT COUNT(*) FROM fg WHERE NOT EXISTS(SELECT 1 FROM fg_reject x WHERE x.c=fg.c AND x.p=fg.p)').fetchone()[0]
    counts['discovery_excluded_unique_positions'] = db.execute('SELECT COUNT(*) FROM fg_reject').fetchone()[0]
    db.commit()
    return dict(counts)


def align(db, out, cfg):
    dcounts = load_discovery(db, out, cfg)
    db.execute("""UPDATE v SET status=CASE WHEN EXISTS(SELECT 1 FROM fg_pos x WHERE x.c=v.c38 AND x.p=v.p38 AND x.row_count>1)
        THEN 'DISCOVERY_POSITION_EXCLUDED' WHEN EXISTS(SELECT 1 FROM fg_reject x WHERE x.c=v.c38 AND x.p=v.p38) THEN 'DISCOVERY_POSITION_EXCLUDED'
        ELSE 'NO_DISCOVERY_POSITION' END WHERE status='READY' AND NOT EXISTS(SELECT 1 FROM fg x WHERE x.c=v.c38 AND x.p=v.p38 AND NOT EXISTS(SELECT 1 FROM fg_reject r WHERE r.c=x.c AND r.p=x.p))""")
    db.execute("""UPDATE v SET status=CASE WHEN v.ea38=f.alt AND v.oa38=f.ref THEN 'ALIGNED_DIRECT'
        WHEN v.ea38=f.ref AND v.oa38=f.alt THEN 'ALIGNED_SWAP' ELSE 'ALLELE_MISMATCH' END
        FROM fg f WHERE v.status='READY' AND v.c38=f.c AND v.p38=f.p AND NOT EXISTS(SELECT 1 FROM fg_reject r WHERE r.c=f.c AND r.p=f.p)""")
    db.commit()
    cols = ['source_row','chromosome','position_grch38','discovery_variant_id','effect_allele','other_allele','bbj_beta_aligned','bbj_se','bbj_p','bbj_eaf_aligned','bbj_n_total','finngen_beta','finngen_se','finngen_p','finngen_eaf','source_position_grch37','source_effect_allele','source_other_allele','lift_strand','alignment','bbj_beta_source','bbj_eaf_source']
    query = '''SELECT v.id,v.c38,v.p38,f.vid,f.alt,f.ref,CASE WHEN v.status='ALIGNED_SWAP' THEN -v.beta ELSE v.beta END,v.se,v.p,
        CASE WHEN v.status='ALIGNED_SWAP' THEN 1-v.eaf ELSE v.eaf END,NULL,f.beta,f.se,f.pval,f.eaf,v.p37,v.ea37,v.oa37,v.strand,v.status,v.beta,v.eaf
        FROM v JOIN fg f ON f.c=v.c38 AND f.p=v.p38 WHERE v.status IN ('ALIGNED_DIRECT','ALIGNED_SWAP') ORDER BY v.c38,v.p38'''
    write_table(out/'bbj_finngen_aiht.aligned.grch38.tsv.gz', cols, db.execute(query))
    write_table(out/'bbj_mapping_alignment_audit.tsv.gz', ['source_row','chromosome_grch37','position_grch37','source_effect_allele','source_other_allele','chromosome_grch38','position_grch38','lift_strand','status'],
                db.execute('SELECT id,c37,p37,ea37,oa37,c38,p38,strand,status FROM v ORDER BY id'))
    counts = dcounts | dict(db.execute('SELECT status,COUNT(*) FROM v GROUP BY status'))
    write_table(out/'alignment_counts.tsv', ['status','records'], sorted(counts.items()))
    return counts


def finish(db, out, cfg, provenance, state):
    checks = []
    def check(name, ok): checks.append({'check': name, 'passed': bool(ok)})
    status_counts = dict(db.execute('SELECT status,COUNT(*) FROM v GROUP BY status'))
    check('official_md5', digest(ROOT/cfg['bbj_input'], 'md5') == cfg['bbj_official_md5'])
    check('raw_row_count', state['clean']['raw_rows'] == cfg['bbj_expected_rows'])
    check('x_row_count', state['clean']['source_x_rows'] == cfg['bbj_expected_x_rows'])
    excluded = sum(value for key, value in state['clean'].items() if key not in {'raw_rows', 'retained', 'candidate_before_source_duplicates', 'source_x_rows'})
    check('raw_exclusion_retained_partition', state['clean']['raw_rows'] == excluded + state['clean']['retained'])
    check('numeric_source_order_checked', (out/'source_order_check.tsv').exists())
    check('source_clean_unique_positions', db.execute("SELECT COUNT(*) FROM v WHERE status NOT LIKE 'MULTIPLE_RECORDS_AT_SOURCE_POSITION'").fetchone()[0] >= state['clean']['retained'])
    check('mapping_partition', sum(state['lift'].values()) == state['clean']['retained'])
    check('final_partition', sum(status_counts.values()) == state['clean']['retained'])
    check('no_pending_states', db.execute("SELECT COUNT(*) FROM v WHERE status IN ('CANDIDATE','READY','ROUNDTRIP_PENDING')").fetchone()[0] == 0)
    check('aligned_nonempty', sum(v for k,v in status_counts.items() if k.startswith('ALIGNED_')) > 0)
    check('aligned_unique_target_positions', db.execute("SELECT COUNT(*) FROM (SELECT c38,p38 FROM v WHERE status LIKE 'ALIGNED_%' GROUP BY c38,p38 HAVING COUNT(*)>1)").fetchone()[0] == 0)
    check('aligned_exact_roundtrip', db.execute("SELECT COUNT(*) FROM v LEFT JOIN rv_one r ON r.id=v.id WHERE v.status LIKE 'ALIGNED_%' AND (r.id IS NULL OR r.c!=CAST(v.c37 AS TEXT) OR r.p!=v.p37 OR r.width!=1 OR r.strand!='+')").fetchone()[0] == 0)
    check('no_hla_in_either_build', db.execute("SELECT COUNT(*) FROM v WHERE status LIKE 'ALIGNED_%' AND ((c37=6 AND p37 BETWEEN 25000000 AND 34000000) OR (c38=6 AND p38 BETWEEN 25000000 AND 34000000))").fetchone()[0] == 0)
    check('no_palindromic_aligned', db.execute("SELECT COUNT(*) FROM v WHERE status LIKE 'ALIGNED_%' AND ea37||oa37 IN ('AT','TA','CG','GC')").fetchone()[0] == 0)
    with gzip.open(out/'bbj_finngen_aiht.aligned.grch38.tsv.gz', 'rt') as handle:
        aligned_rows = 0
        n_all_na = True
        for row in csv.DictReader(handle, delimiter='\t'):
            aligned_rows += 1
            n_all_na = n_all_na and row['bbj_n_total'] == 'NA'
    check('bbj_n_total_all_na', n_all_na)
    check('aligned_output_row_count', aligned_rows == db.execute("SELECT COUNT(*) FROM v WHERE status LIKE 'ALIGNED_%'").fetchone()[0])
    for path, expected in provenance.items(): check('input_unchanged:' + path, digest(ROOT/path) == expected)
    (out/'validation.json').write_text(json.dumps(checks, indent=2) + '\n')
    decision = {'completed_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS_HARMONIZATION_HOLD_DOWNSTREAM_ANALYSES' if all(x['passed'] for x in checks) else 'HOLD_VALIDATION',
        'source_records': state['clean']['raw_rows'], 'clean_records': state['clean']['retained'], 'aligned_records': sum(v for k,v in status_counts.items() if k.startswith('ALIGNED_')),
        'alignment_statuses': status_counts, 'checks_passed': sum(x['passed'] for x in checks), 'checks_total': len(checks),
        'liftover_executed': True, 'roundtrip_verified': True, 'replication_executed': False, 'heterogeneity_executed': False, 'coloc_executed': False,
        'metadata_only': {'bbj_cohort_n_cases': cfg['bbj_cohort_n_cases_metadata'], 'bbj_cohort_n_controls': cfg['bbj_cohort_n_controls_metadata']},
        'holds': ['No per-variant BBJ N, cases, controls, MAC, or Rsq supplied', 'No replication, heterogeneity, LD, or colocalisation analyses', 'EA/OA are not asserted to be BBJ reference/alternate alleles', 'Target matching requires exact FinnGen ALT/REF concordance after validated lift strand handling']}
    (out/'decision.json').write_text(json.dumps(decision, indent=2) + '\n')
    assert all(x['passed'] for x in checks), json.dumps(checks, indent=2)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', default='outcome_data/t16/harmonization_20260904')
    args = parser.parse_args(); out = (ROOT / args.output).resolve()
    assert out.is_relative_to(ROOT/'outcome_data/t16')
    if out.exists(): raise FileExistsError(f'Refusing to overwrite existing output: {out}')
    out.mkdir(parents=True)
    cfgpath = ROOT/'config/t16_harmonization.json'; cfg = json.loads(cfgpath.read_text())
    inputs = [cfgpath, Path(__file__).resolve(), ROOT/cfg['bbj_input'], ROOT/cfg['discovery_input'], ROOT/cfg['forward_chain'], ROOT/cfg['reverse_chain'], ROOT/cfg['liftover_executable'], ROOT/cfg['t17_reference_script']]
    provenance = {str(path.relative_to(ROOT)): digest(path) for path in inputs}
    state = {'provenance': provenance}
    (out/'checkpoint.json').write_text(json.dumps(state, indent=2) + '\n')
    db = sqlite3.connect(out/'harmonization.sqlite'); db.execute('PRAGMA cache_size=-262144'); db.execute('PRAGMA journal_mode=WAL'); db.execute('PRAGMA synchronous=NORMAL')
    state['clean'] = prepare_source(db, out, cfg)
    # Store verified preflight values explicitly instead of deriving from an output row count.
    state['lift'] = lift(db, out, cfg); state['align'] = align(db, out, cfg)
    (out/'checkpoint.json').write_text(json.dumps(state, indent=2) + '\n')
    finish(db, out, cfg, provenance, state)
    write_table(out/'provenance.tsv', ['path','sha256'], provenance.items())
    db.execute('PRAGMA wal_checkpoint(TRUNCATE)'); db.close()


if __name__ == '__main__': main()
