#!/usr/bin/env python3
"""Independent small-artifact validation of T18, without importing its producer."""
import argparse
import collections
import csv
import hashlib
import json
import math
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output')
    args = ap.parse_args()
    out = Path(args.output).resolve()
    assert not (out / 'independent_synthesis_validation.json').exists()
    checks = []

    def check(name, condition):
        checks.append({'check': name, 'passed': bool(condition)})

    manifest = json.loads((out / 'input_manifest.json').read_text())
    check('all_current_boundary_hashes_match', all(sha(ROOT / p) == h for p, h in manifest['sha256'].items()))
    source = ROOT / 'outcome_data/external_validation/run_20260904/priority_variant_validation.tsv'
    check('source_statistics_byte_identical', sha(source) == sha(out / source.name))
    raw = rows(source)
    lookup = {(r['signal_unit_id'], r['trait']): r for r in raw}
    summary = rows(out / 'class_summary.tsv')
    check('all_five_classes_by_two_traits', len(summary) == 10 and len({(r['discovery_class'], r['trait']) for r in summary}) == 10)
    for r in summary:
        sub = [v for v in raw if v['five_class_discovery'] == r['discovery_class'] and v['trait'] == r['trait']]
        check('counts_' + r['discovery_class'] + '_' + r['trait'],
              int(r['planned_variant_tests']) == len(sub) and
              int(r['testable_variant_tests']) == sum(v['availability'] == 'aligned_exact_target' for v in sub) and
              int(r['same_direction_association_support']) == sum(v['classification'] == 'same_direction_association_support' for v in sub) and
              int(r['opposite_direction_association_support']) == sum(v['classification'] == 'opposite_direction_association_support' for v in sub) and
              int(r['inconclusive_low_precision']) == sum(v['classification'] == 'inconclusive_low_precision_OR1_10_scenario' for v in sub) and
              int(r['untestable']) == sum(v['availability'] != 'aligned_exact_target' for v in sub))
    check('no_groupwise_denominator_reduction', sum(int(r['planned_variant_tests']) for r in summary) == 20)
    units = rows(ROOT / 'outcome_data/t07/genetic_signal_five_classification.tsv')
    count = collections.Counter(r['five_class_primary'] for r in units)
    check('discovery_class_counts', all(int(r['discovery_units_total']) == count[r['discovery_class']] for r in summary))
    plot = rows(out / 'forest_data.tsv')
    keys = {(r['signal_unit_id'], r['trait'], r['cohort']) for r in plot}
    check('plot_full_grid', len(plot) == len(keys) == 40 and keys == {(sid, tr, c) for sid, tr in lookup for c in ['FinnGen', 'BBJ']})
    for r in plot:
        source_row = lookup[(r['signal_unit_id'], r['trait'])]
        prefix = 'external' if r['cohort'] == 'BBJ' else 'discovery'
        b, s = source_row[prefix + '_beta'], source_row[prefix + '_se']
        if b == 'NA':
            valid = r['available'] == 'False' and all(r[k] == 'NA' for k in ['beta', 'se', 'odds_ratio', 'lower95', 'upper95', 'source_p'])
        else:
            b, s = float(b), float(s)
            expected = {'beta': b, 'se': s, 'odds_ratio': math.exp(b),
                        'lower95': math.exp(b - 1.959963984540054 * s),
                        'upper95': math.exp(b + 1.959963984540054 * s)}
            valid = r['available'] == 'True' and all(math.isclose(float(r[k]), v, rel_tol=1e-12, abs_tol=1e-12) for k, v in expected.items())
            valid = valid and r['source_p'] == source_row['external_source_p' if prefix == 'external' else 'discovery_p']
        valid = valid and all(r[k] == source_row[k] for k in ['variant_id', 'effect_allele', 'other_allele', 'nearest_genes'])
        valid = valid and r['external_classification'] == source_row['classification']
        check('plot_' + '_'.join([r['signal_unit_id'], r['trait'], r['cohort']]), valid)
    check('missing_not_zero', sum(r['available'] == 'False' for r in plot) == 10)
    pair = rows(out / 'representative_evidence.tsv')
    check('pair_ledger_count', len(pair) == 10 and len({r['signal_unit_id'] for r in pair}) == 10)
    pair_ok = True
    for r in pair:
        for p, trait in [('graves', 'Graves'), ('aiht', 'AIHT')]:
            for field in ['classification', 'external_beta', 'external_source_p', 'direction', 'external_alignment', 'external_lift_strand', 'heterogeneity_p_bonferroni']:
                pair_ok &= r[p + '_' + field] == lookup[(r['signal_unit_id'], trait)][field]
    check('pair_ledger_exact_fields', pair_ok)
    opposite = rows(out / 'shared_opposite_audit.tsv')
    check('opposite_both_targets_and_holds', {r['signal_unit_id'] for r in opposite} == {'GSU056', 'GSU254'} and
          all(r['different_signal_excluded'] == r['ancestry_ld_difference_excluded'] == r['cross_disease_covariance_available'] == 'False' for r in opposite))
    for trait in ['Graves', 'AIHT']:
        png = (out / ('forest_' + trait + '.png')).read_bytes()
        check('png_' + trait, png[:8] == b'\x89PNG\r\n\x1a\n' and struct.unpack('>II', png[16:24]) == (2600, 1600))
        svg = ET.parse(out / ('forest_' + trait + '.svg')).getroot()
        check('svg_' + trait, svg.tag.endswith('svg') and len(list(svg.iter())) > 100)
    numeric = rows(out / 'independent_statistics_validation.tsv')
    check('fresh_R_29_pass', len(numeric) == 29 and all(r['passed'] == 'TRUE' for r in numeric))
    decision = json.loads((out / 'decision.json').read_text())
    check('decision_scope', decision['variant_level_synthesis_completed'] and decision['planned_tests'] == 20 and decision['testable_tests'] == 10 and
          not any(decision[k] for k in ['full_t18_complete', 'signal_level_validation_completed', 'cross_disease_difference_executed', 'manuscript_changed', 't19_started']))
    ledger = rows(out / 'claim_disposition.tsv')
    check('claim_ledger_scope', len(ledger) == 10 and any(r['action'] == 'ALLOW_VARIANT_ASSOCIATION_ONLY' for r in ledger) and any(r['action'] == 'DO_NOT_CLAIM' for r in ledger))
    result = {'status': 'PASS' if all(c['passed'] for c in checks) else 'FAIL',
              'checks_passed': sum(c['passed'] for c in checks), 'checks_total': len(checks),
              'validator_sha256': sha(Path(__file__)), 'checks': checks,
              'output_sha256': {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
    (out / 'independent_synthesis_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'checks_passed', 'checks_total']}))
    assert result['status'] == 'PASS'


if __name__ == '__main__':
    main()
