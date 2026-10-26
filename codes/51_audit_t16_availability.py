#!/usr/bin/env python3
"""Rebuild T16 availability evidence from frozen public snapshots; no network."""
import csv
import hashlib
import json
import zlib
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outcome_data/t16'


class PortalRows(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell, self.links = [], None, None, []

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row, self.links = [], []
        elif tag == 'td' and self.row is not None:
            self.cell = []
        elif tag == 'a' and self.row is not None:
            self.links.append(dict(attrs).get('href', ''))

    def handle_data(self, text):
        if self.cell is not None:
            self.cell.append(text)

    def handle_endtag(self, tag):
        if tag == 'td' and self.cell is not None:
            self.row.append(''.join(self.cell).strip())
            self.cell = None
        elif tag == 'tr' and self.row is not None:
            for link in self.links:
                if '/dataset_hypothyroidism/' in link:
                    self.rows.append({'trait': self.row[1], 'stratum': self.row[2],
                                      'size_mb_as_reported': self.row[3],
                                      'url': urljoin('https://genepi.med.uni-greifswald.de', link),
                                      'non_ukb_split': False})
            self.row = None


def write(name, rows):
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)


def main():
    config = json.loads((ROOT / 'config/t16_availability.json').read_text())
    portal = (OUT / 'sources/thyroidomics_datasets.html').read_text()
    p = PortalRows(); p.feed(portal)
    assert len(p.rows) == 5, 'Portal structure or public release changed; review before regenerating'
    write('t16_thyroidomics_public_files.tsv', p.rows)
    write('t16_availability_registry.tsv', config['candidates'])
    sources = config['sources']
    write('t16_sources.tsv', sources)
    checks = [
        {'check': 'five_public_thyroidomics_strata', 'passed': len(p.rows) == 5},
        {'check': 'no_observed_non_ukb_split', 'passed': not any('ukb' in r['url'].lower() for r in p.rows)},
        {'check': 'no_eligible_dataset_without_file', 'passed': all(
            not r['ready_for_replication'] for r in config['candidates'])},
        {'check': 'no_external_request_sent', 'passed': config['external_request_sent'] is False},
        {'check': 'frozen_exclusions_retained', 'passed': config['excluded_cohorts'] == ['FINNGEN', 'UKB']},
        {'check': 'mvp_example_not_aiht', 'passed': 'Phe_244.AMR.GIA' in
            (OUT / 'sources/mvp_pha6812.html').read_text()},
    ]
    bbj = json.loads((OUT / 'sources/GCST90018635.json').read_text())
    metadata = (OUT / 'sources/GCST90018635_meta.yaml').read_text()
    prefix = (OUT / 'sources/GCST90018635_prefix.gz.part').read_bytes()
    decompressed = zlib.decompressobj(31).decompress(prefix).decode()
    header = decompressed.splitlines()[0].split('\t')
    checks += [
        {'check': 'bbj_only_hashimoto_accession', 'passed': bbj['cohort'] == 'BBJ' and
         bbj['diseaseTrait']['trait'] == 'Hashimoto thyroiditis'},
        {'check': 'bbj_sample_counts', 'passed': bbj['initialSampleSize'] ==
         '537 East Asian ancestry cases, 172,656 East Asian ancestry controls'},
        {'check': 'bbj_grch37_metadata', 'passed': 'genome_assembly: GRCh37' in metadata},
        {'check': 'bbj_actual_header', 'passed': {'chromosome', 'base_pair_location',
         'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value'} <= set(header)},
        {'check': 'only_prefix_fetched', 'passed': len(prefix) == 65536 and
         'Content-Range: bytes 0-65535/526237307' in (OUT / 'sources/GCST90018635_range_headers.txt').read_text()},
        {'check': 'upstream_md5_consistent', 'passed': '91fc15fd74923c4f01ce31f2c87cb2ef' in metadata and
         '91fc15fd74923c4f01ce31f2c87cb2ef' in (OUT / 'sources/GCST90018635_md5sum.txt').read_text()},
    ]
    (OUT / 't16_bbj_header_preview.txt').write_text('\n'.join(decompressed.splitlines()[:6]) + '\n')
    addendum = []
    with (ROOT / 'outcome_data/t15/t15_external_cohort_registry.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            cohorts = set(row['cohorts'].split(';'))
            relation = ('UNCONFIRMED' if row['cohort_membership_complete'] != 'True' else
                        'PARTIAL_OVERLAP' if 'BBJ' in cohorts else 'NO_SHARED_COHORT_IDENTIFIED')
            # Unknown cross-enrolment in a consumer cohort remains unconfirmed.
            if '23ANDME' in cohorts and 'BBJ' not in cohorts:
                relation = 'UNCONFIRMED'
            addendum.append({'dataset_a': 'BBJ_2021_HASHIMOTO_GCST90018635',
                             'dataset_b': row['dataset_id'], 'overlap_class': relation,
                             'shared_cases_n': None, 'shared_controls_n': None,
                             'basis': 'Cohort sampling frame only, not participant-level linkage'})
    write('t16_bbj_overlap_addendum.tsv', addendum)
    write('t16_validation.tsv', checks)
    provenance = []
    paths = [ROOT / 'config/t16_availability.json', Path(__file__).resolve(),
             ROOT / 'outcome_data/t15/t15_external_cohort_registry.tsv']
    paths += [ROOT / s['snapshot'] for s in sources if s['snapshot']]
    for path in dict.fromkeys(paths):
        provenance.append({'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size,
                           'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    write('t16_provenance.tsv', provenance)
    decision = {'audit_date': config['audit_date'], 'status': 'PASS_AVAILABILITY_BBJ_HOLD_FULL_QC_AND_EUR',
                'scope': 'Public availability audit only; bounded candidates, not exhaustive',
                'candidate_records': len(config['candidates']),
                'checks_passed': sum(x['passed'] for x in checks), 'checks_total': len(checks),
                'replication_executed': False, 'external_request_sent': False,
                'ready_for_acquisition': ['BBJ_2021_HASHIMOTO_GCST90018635'],
                'full_gwas_downloaded': False,
                'note': 'Public independent BBJ Hashimoto file confirmed; only prefix inspected. EUR input and full QC remain held.'}
    (OUT / 't16_decision.json').write_text(json.dumps(decision, indent=2) + '\n')
    print(json.dumps(decision, indent=2))
    assert all(x['passed'] for x in checks)


if __name__ == '__main__':
    main()
