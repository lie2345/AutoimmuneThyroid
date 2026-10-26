#!/usr/bin/env python3
"""Offline acquisition validation and descriptive QC, never replication analysis.

Preserves the official archive; extracts only the combined-sex autosome/X members.
No liftover, effect harmonisation, locus testing, or LD substitution is performed.
"""
import collections
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path
import shutil
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'raw_data/external_gwas/bbj_2019_graves'
OUT = ROOT / 'outcome_data/t17'
ARCHIVE = RAW / 'hum0014.v17.GD.v1.zip'
EXPECTED_BYTES = 2245854863
MEMBERS = ['GD.auto.rsq07.mac10.txt.gz', 'GD.chrx.rsq07.mac10.txt.gz']
NS = 212453


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def write_tsv(name, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def decode_opengwas():
    j = json.loads((OUT / 'sources/opengwas_gwasinfo.json').read_text())
    fields = ['trait', 'build', 'category', 'subcategory', 'population', 'sex',
              'author', 'year', 'ontology', 'unit', 'sample_size', 'consortium', 'mr', 'priority']
    a, b = j['datasets_compressed']['bbj-a-123']
    d = {}
    for k, v in zip(fields, a):
        coding = j['majority_fields_and_coding'].get(k, [])
        d[k] = coding[v] if isinstance(v, int) and 0 <= v < len(coding) else v
    return dict(d, **b)


def sample_counts():
    """Read publisher cells without modifying the workbook or its extensions."""
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(OUT / 'sources/hum0014_sample_size.xlsx') as z:
        strings = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        cells = {}
        for c in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//m:c', ns):
            v = c.find('m:v', ns)
            if v is not None:
                cells[c.attrib['r']] = strings[int(v.text)] if c.get('t') == 's' else v.text
    assert cells['A24'] == "Graves' disease" and cells['C24'] == 'hum0014.v17.GD.v1'
    return {k: int(cells[k]) for k in ['D24', 'E24', 'F24', 'G24', 'H24', 'I24']}


def row_errors(r, x=False):
    errors = []
    def number(k, predicate):
        try:
            v = float(r[k])
            if not math.isfinite(v) or not predicate(v):
                errors.append(k)
        except (ValueError, KeyError):
            errors.append(k)
    number('POS', lambda v: v > 0 and v == int(v))
    number('BETA', lambda v: True)
    number('SE', lambda v: v > 0)
    number('p.value', lambda v: 0 <= v <= 1)
    a, b = r.get('Allele1', ''), r.get('Allele2', '')
    def valid_allele(v):
        return bool(v) and (set(v) <= set('ACGT') or re.fullmatch(r'<[A-Z0-9:]+>', v))
    if not valid_allele(a) or not valid_allele(b) or a == b:
        errors.append('alleles')
    if r.get('SNPID') != '_'.join(r.get(k, '') for k in ['CHR', 'POS', 'Allele1', 'Allele2']):
        errors.append('SNPID')
    if x:
        for suffix in ['x', 'y']:
            for base in ['AF.Cases', 'AF.Controls']:
                number(base + '.' + suffix, lambda v: 0 <= v <= 1)
            number('MAC.' + suffix, lambda v: v >= 0)
            number('N.' + suffix, lambda v: v > 0)
        for k in ['Rsq.male', 'Rsq.female']:
            number(k, lambda v: 0 <= v <= 1)
        try:
            if float(r['N.x']) + float(r['N.y']) != NS:
                errors.append('N.total')
        except (KeyError, ValueError):
            errors.append('N.total')
    else:
        # Variant-specific N is legitimate. Never overwrite it with the cohort N.
        number('N', lambda v: 0 < v <= NS and v == int(v))
        number('Rsq', lambda v: 0.7 <= v <= 1)
        number('MAC', lambda v: v >= 10)
        number('AC_Allele2', lambda v: 0 <= v <= 2 * NS)
        for k in ['AF_Allele2', 'AF.Cases', 'AF.Controls']:
            number(k, lambda v: 0 <= v <= 1)
        if r.get('Is.SPA.converge') not in ['0', '1']:
            errors.append('SPA_flag')
    return errors


def scan(path):
    x = '.chrx.' in path.name
    counts = collections.Counter()
    chr_counts = collections.Counter()
    n_counts = collections.Counter()
    duplicate_ids = collections.Counter()
    errors = collections.Counter()
    examples = []
    # SNPID uniqueness is checked within contiguous chromosome groups. A return
    # to an earlier chromosome is a hard error, not silently dropped coverage.
    seen, finished, previous_chr, previous_pos = set(), set(), None, None
    with gzip.open(path, 'rt') as f:
        header = next(f).split()
        required = {'CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'BETA', 'SE', 'p.value'}
        assert required <= set(header), (path.name, header)
        for i, line in enumerate(f, 1):
            counts['rows'] += 1
            fields = line.split()
            if len(fields) != len(header):
                errors['field_count'] += 1
                continue
            r = dict(zip(header, fields))
            e = row_errors(r, x=x)
            for k in e:
                errors[k] += 1
            if e and len(examples) < 10:
                examples.append({'row': i, 'id': r.get('SNPID'), 'errors': e})
            c = r['CHR']
            chr_counts[c] += 1
            if c != previous_chr:
                if previous_chr is not None:
                    finished.add(previous_chr)
                if c in finished:
                    errors['noncontiguous_chromosome'] += 1
                seen.clear()
                previous_pos = None
                previous_chr = c
            if r['SNPID'] in seen:
                counts['duplicate_variant_id'] += 1
                duplicate_ids[r['SNPID']] += 1
            seen.add(r['SNPID'])
            try:
                pos = int(r['POS'])
                if previous_pos is not None and pos < previous_pos:
                    counts['position_decreases'] += 1
                previous_pos = pos
                if c == '6' and 25000000 <= pos <= 34000000:
                    counts['hla_25_34mb_grch37'] += 1
            except ValueError:
                pass
            a, b = r['Allele1'], r['Allele2']
            if not set(a + b) <= set('ACGT'):
                counts['symbolic_alleles_requires_separate_handling'] += 1
            counts['snv' if len(a) == len(b) == 1 else 'indel_or_complex'] += 1
            if a + b in ['AT', 'TA', 'CG', 'GC']:
                counts['palindromic_snv'] += 1
            if r['p.value'] in ['0', '0.0']:
                counts['p_zero_underflow'] += 1
            if not x and r.get('Is.SPA.converge') == '0':
                counts['spa_not_converged'] += 1
            if not x:
                n_counts[r['N']] += 1
                if r['N'] != str(NS):
                    counts['variant_specific_n_differs_from_cohort'] += 1
            if x:
                if any(float(r[k]) < 0.7 for k in ['Rsq.male', 'Rsq.female']):
                    counts['x_any_sex_rsq_below_0_7'] += 1
                if any(float(r[k]) < 10 for k in ['MAC.x', 'MAC.y']):
                    counts['x_any_sex_mac_below_10'] += 1
            if i % 2000000 == 0:
                print(f'{path.name}: {i:,} rows scanned', flush=True)
    return {'file': path.name, 'header': header, 'counts': dict(counts),
            'chromosomes': dict(chr_counts), 'n_counts': dict(n_counts),
            'duplicate_ids_extra_records': dict(duplicate_ids),
            'errors': dict(errors), 'examples': examples}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    candidate = ARCHIVE if ARCHIVE.exists() else ARCHIVE.with_suffix('.zip.part')
    if not candidate.exists() or candidate.stat().st_size != EXPECTED_BYTES:
        raise SystemExit('HOLD: complete official ZIP not yet available; no QC PASS written')
    print('Checking CRC of every archive member...', flush=True)
    with zipfile.ZipFile(candidate) as z:
        assert z.testzip() is None, 'ZIP CRC failed'
        inventory = [{'member': i.filename, 'bytes': i.file_size,
                      'compressed_bytes': i.compress_size, 'crc32': f'{i.CRC:08x}',
                      'selected': Path(i.filename).name in MEMBERS} for i in z.infolist()]
        for name in MEMBERS:
            member = 'hum0014.v17.GD.v1/' + name
            target = RAW / name
            # Existing data are reused only after checking against the archive.
            if target.exists():
                import zlib
                crc = 0
                with target.open('rb') as f:
                    for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
                        crc = zlib.crc32(b, crc)
                assert target.stat().st_size == z.getinfo(member).file_size
                assert crc == z.getinfo(member).CRC
            else:
                with z.open(member) as src, target.with_suffix('.gz.part').open('wb') as dst:
                    shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
                target.with_suffix('.gz.part').rename(target)
    if candidate != ARCHIVE:
        candidate.rename(ARCHIVE)
    write_tsv('t17_archive_members.tsv', inventory)
    results = [scan(RAW / name) for name in MEMBERS]
    meta = decode_opengwas()
    (OUT / 't17_opengwas_metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
    (OUT / 't17_file_qc.json').write_text(json.dumps(results, indent=2) + '\n')
    checks = []
    def check(k, ok, detail):
        checks.append({'check': k, 'status': 'PASS' if ok else 'FAIL', 'detail': detail})
    check('official_archive_size_crc', True, str(EXPECTED_BYTES))
    check('opengwas_metadata', meta['ncase'] == 2176 and meta['ncontrol'] == 210277
          and meta['build'] == 'HG19/GRCh37', '2176/210277; GRCh37')
    samples = sample_counts()
    check('official_sample_table', samples['D24'] == 2176 and samples['G24'] == 210277
          and samples['E24'] + samples['F24'] == samples['D24']
          and samples['H24'] + samples['I24'] == samples['G24'],
          'Sample size!D24:I24; ' + json.dumps(samples))
    for r in results:
        check(r['file'] + '_parse_and_fields', not r['errors'], json.dumps(r['errors']))
        expected_chr = {'X'} if '.chrx.' in r['file'] else {str(i) for i in range(1, 23)}
        check(r['file'] + '_chromosomes', set(r['chromosomes']) == expected_chr,
              ';'.join(r['chromosomes']))
        check(r['file'] + '_nonempty', r['counts']['rows'] > 100000, str(r['counts']['rows']))
    dictionary = (OUT / 'sources/hum0014_README.txt').read_text()
    check('effect_allele_dictionary', 'Effect allele = alternative allele' in dictionary
          and 'Allele2: alternative allele' in dictionary, 'BBJ Allele2, not Allele1')
    check('discovery_build_unchanged', 'genome_build: "GRCh38"' in (ROOT / 'config/analysis.yml').read_text(),
          'GRCh38 discovery must be harmonised to GRCh37 separately')
    total = sum(r['counts']['rows'] for r in results)
    check('source_total_matches_opengwas_count', total == meta['nsnp'],
          f'official rows={total}; OpenGWAS nsnp={meta["nsnp"]}; count agreement is not byte identity')
    write_tsv('t17_validation.tsv', checks)
    prep_holds = [
        {'item': 'duplicate_ID_resolution', 'status': 'HOLD',
         'evidence': json.dumps([r['duplicate_ids_extra_records'] for r in results])},
        {'item': 'symbolic_allele_handling', 'status': 'HOLD',
         'evidence': str(sum(r['counts'].get('symbolic_alleles_requires_separate_handling', 0) for r in results))},
        {'item': 'per_variant_N_preservation', 'status': 'REQUIRED',
         'evidence': json.dumps(results[0]['n_counts']) + '; per-variant case/control split unknown'},
        {'item': 'numeric_sort_before_indexing', 'status': 'REQUIRED',
         'evidence': 'autosomes lexicographic chromosome order; X has position decreases'},
    ]
    write_tsv('t17_variant_preparation_gates.tsv', prep_holds)
    decision = {
        'completed_at_utc': datetime.now(timezone.utc).isoformat(),
        'stage': 'T17 acquisition and preflight only',
        'status': 'PASS_WITH_HOLD' if all(c['status'] == 'PASS' for c in checks) else 'HOLD_QC_REVIEW',
        'checks_passed': sum(c['status'] == 'PASS' for c in checks), 'checks_total': len(checks),
        'source': 'NBDC hum0014.v17.GD.v1; not an OpenGWAS VCF download',
        'effect_allele': 'Allele2', 'genome_build': 'GRCh37',
        'full_archive_downloaded': True, 'combined_sex_files_scanned': True,
        'replication_executed': False, 'liftover_executed': False,
        'holds': ['Duplicate IDs and symbolic alleles require explicit handling before analysis',
                  'Preserve per-variant N; case/control split for N=212451 is unknown',
                  'Discovery GRCh38 vs BBJ GRCh37 harmonisation not executed',
                  'No validated East Asian/Japanese LD panel configured',
                  'Exact clinical phenotype equivalence to FinnGen strict Graves not established',
                  'Independent European Graves source not acquired',
                  'Variant-level filter/sort/normalisation and effect alignment not executed'],
    }
    (OUT / 't17_decision.json').write_text(json.dumps(decision, indent=2) + '\n')
    files = [ARCHIVE, *[RAW / n for n in MEMBERS], Path(__file__).resolve(),
             ROOT / 'config/analysis.yml', ROOT / 'config/traits.tsv',
             ROOT / 'outcome_data/t15/t15_decision.tsv',
             OUT / 'sources/hum0014_README.txt', OUT / 'sources/hum0014_sample_size.xlsx',
             OUT / 'sources/hum0014_v17_GD_v1.html', OUT / 'sources/opengwas_gwasinfo.json']
    write_tsv('t17_provenance.tsv', [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size,
                                    'sha256': sha256(p)} for p in files])
    print(json.dumps(decision, indent=2), flush=True)


if __name__ == '__main__':
    main()
