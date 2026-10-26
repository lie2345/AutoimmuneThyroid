#!/usr/bin/env python3
"""Summarize the verified variant-level boundary; never reclassify T07 signals."""
import argparse
import collections
import csv
import hashlib
import json
import math
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'outcome_data/external_validation/run_20260904'
PLAN = ROOT / 'outcome_data/external_validation/execution_plan_20260904'
CLASSES = ['shared-concordant', 'shared-opposite', 'Graves-predominant', 'AIHT-predominant', 'unresolved']
LABELS = {'same_direction_association_support': '同向变异关联支持',
          'opposite_direction_association_support': '反向变异关联支持',
          'no_primary_support_small_effect_equivalence': '未达关联阈值／小效应等效',
          'inconclusive_low_precision_OR1_10_scenario': '未定／OR1.10情景低精度',
          'inconclusive': '未定', 'untestable_missing_aligned_target': '不可检验'}


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def table(path, rows):
    assert rows
    with path.open('x') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def summarize(rows, targets, universe):
    result = []
    counts = collections.Counter(r['five_class_primary'] for r in universe)
    for cls in CLASSES:
        planned = sum(r['five_class_primary'] == cls for r in targets)
        for trait in ['Graves', 'AIHT']:
            sub = [r for r in rows if r['five_class_discovery'] == cls and r['trait'] == trait]
            tally = collections.Counter(r['classification'] for r in sub)
            result.append(dict(discovery_class=cls, trait=trait, discovery_units_total=counts[cls],
                planned_variant_tests=planned, testable_variant_tests=sum(r['availability'] == 'aligned_exact_target' for r in sub),
                same_direction_association_support=tally['same_direction_association_support'],
                opposite_direction_association_support=tally['opposite_direction_association_support'],
                inconclusive_low_precision=tally['inconclusive_low_precision_OR1_10_scenario'],
                other_testable_classifications=sum(tally[k] for k in ['inconclusive', 'no_primary_support_small_effect_equivalence']),
                untestable=tally['untestable_missing_aligned_target'],
                interpretation='NOT_IN_THIS_VALIDATION_SCOPE' if not planned else 'VARIANT_ONLY_SIGNAL_CLASS_NOT_REVALIDATED'))
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output).resolve()
    assert out.is_relative_to(ROOT / 'outcome_data/t18'), 'Output must be a new T18 child directory'
    assert not out.exists(), 'Refusing to overwrite existing output'
    inputs = [SOURCE / x for x in ['priority_variant_validation.tsv', 'decision.json', 'provenance.tsv', 'independent_statistics_validation.tsv']]
    inputs += [PLAN / x for x in ['plan.json', 'targets.tsv', 'freeze.json', 'amendment.json']]
    inputs += [ROOT / x for x in ['outcome_data/t07/genetic_signal_five_classification.tsv',
        'outcome_data/t07/priority_signal_units.tsv', 'codes/59_verify_external_statistics.R',
        'codes/63_synthesize_t18.py', 'codes/64_plot_t18_forest.R']]
    before = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    frozen = json.loads((PLAN / 'freeze.json').read_text())['sha256']
    # Verify frozen small contracts, not multi-GB upstream data already fully audited.
    for p in [PLAN / 'plan.json', PLAN / 'targets.tsv', PLAN / 'amendment.json',
              ROOT / 'outcome_data/t07/priority_signal_units.tsv', ROOT / 'codes/57_run_external_validation.py']:
        assert sha(p) == frozen[str(p.relative_to(ROOT))], f'Frozen contract changed: {p}'
    rows, targets = read(inputs[0]), read(PLAN / 'targets.tsv')
    universe = read(ROOT / 'outcome_data/t07/genetic_signal_five_classification.tsv')
    original_decision = json.loads((SOURCE / 'decision.json').read_text())
    assert original_decision['status'] == 'PASS_FEASIBLE_VARIANT_VALIDATION_HOLD_SIGNAL_LEVEL_AND_T18_SYNTHESIS'
    assert all(r['passed'] == 'TRUE' for r in read(SOURCE / 'independent_statistics_validation.tsv'))
    assert len(rows) == 20 and len(targets) == 10 and len(universe) == 304
    targetmap = {r['signal_unit_id']: r for r in targets}
    unitmap = {r['signal_unit_id']: r for r in universe}
    assert len(targetmap) == 10 and len(unitmap) == 304
    assert {(r['signal_unit_id'], r['trait']) for r in rows} == {(s, t) for s in targetmap for t in ['Graves', 'AIHT']}
    for r in rows:
        t, u = targetmap[r['signal_unit_id']], unitmap[r['signal_unit_id']]
        assert r['five_class_discovery'] == t['five_class_primary'] == u['five_class_primary']
        assert r['variant_id'] == t['representative_variant'] == u['representative_variant']
        assert (r['effect_allele'], r['other_allele']) == (t['alt'], t['ref']) == (u['alt'], u['ref'])
        assert r['classification'] in LABELS
        assert r['discovery_class_external_status'] == 'NOT_REVALIDATED_NO_EXTERNAL_SIGNAL_LEVEL_EVIDENCE'
    actual_counts = collections.Counter((r['trait'], r['classification']) for r in rows)
    assert actual_counts == collections.Counter({(r['trait'], r['classification']): r['n'] for r in original_decision['counts']})
    assert {(r['signal_unit_id'], r['trait']) for r in rows if r['classification'] == 'same_direction_association_support'} == {('GSU052', 'Graves'), ('GSU073', 'Graves'), ('GSU231', 'Graves')}
    available = [r for r in rows if r['availability'] == 'aligned_exact_target']
    assert all(float(r['heterogeneity_p_bonferroni']) > .05 and r['equivalence_flag'] == 'False' for r in available)
    assert all(float(r['precision_power_or1_10']) < .8 for r in available)
    out.mkdir(parents=True)
    dump(out / 'input_manifest.json', {'created_at_utc': datetime.now(timezone.utc).isoformat(),
         'sha256': before, 'boundary': 'Verified T16/T17 result boundary; full raw/database audit not repeated in T18'})
    # Byte-identical snapshot permits independent R recomputation without modifying prior results.
    shutil.copyfile(SOURCE / 'priority_variant_validation.tsv', out / 'priority_variant_validation.tsv')
    with (out / 'independent_statistics.log').open('x') as log:
        subprocess.run(['Rscript', str(ROOT / 'codes/59_verify_external_statistics.R'), str(out)], stdout=log, stderr=subprocess.STDOUT, check=True)
    checks = read(out / 'independent_statistics_validation.tsv')
    assert len(checks) == 29 and all(r['passed'] == 'TRUE' for r in checks)
    summary = summarize(rows, targets, universe)
    table(out / 'class_summary.tsv', summary)
    # Full representative-level pair ledger preserves both phenotypes without a cross-disease test.
    pairs, reverse = [], []
    by = {(r['signal_unit_id'], r['trait']): r for r in rows}
    for sid, t in targetmap.items():
        g, a = by[(sid, 'Graves')], by[(sid, 'AIHT')]
        pair = dict(signal_unit_id=sid, discovery_class=t['five_class_primary'], variant_id=t['representative_variant'],
                    nearest_genes=t['representative_nearest_genes'], effect_allele=t['alt'], other_allele=t['ref'])
        for name, r in [('graves', g), ('aiht', a)]:
            for key in ['classification', 'external_beta', 'external_source_p', 'direction', 'external_alignment', 'external_lift_strand', 'heterogeneity_p_bonferroni']:
                pair[name + '_' + key] = r[key]
        pair['external_signal_class_status'] = 'HOLD_NO_EAS_LD_OR_SIGNAL_LEVEL_EVIDENCE'
        pairs.append(pair)
        if t['five_class_primary'] == 'shared-opposite':
            reverse.append({**pair,
                'allele_orientation_evidence': 'PASS_FULL_UPSTREAM_HARMONIZATION_AUDIT' if g['availability'] == a['availability'] == 'aligned_exact_target' else 'UNTESTABLE_PRIMARY_VARIANT_MISSING',
                'external_two_trait_sign_pattern': 'same_sign_descriptive_only' if g['availability'] == a['availability'] == 'aligned_exact_target' and float(g['external_beta']) * float(a['external_beta']) > 0 else 'unavailable_or_other',
                'different_signal_excluded': False, 'ancestry_ld_difference_excluded': False,
                'cross_disease_covariance_available': False, 'interpretation': 'NO_EXTERNAL_SHARED_OPPOSITE_VALIDATION_NO_REFUTATION'})
    table(out / 'representative_evidence.tsv', pairs)
    table(out / 'shared_opposite_audit.tsv', reverse)
    claims = [
        dict(claim='T07 discovery five-class results', action='RETAIN_DISCOVERY_SCOPE', evidence='304 discovery units; original T07 outputs unchanged', limitation='None of the signal classes is externally revalidated'),
        dict(claim='Graves associations at GSU052, GSU073, GSU231', action='ALLOW_VARIANT_ASSOCIATION_ONLY', evidence='Same-trait same-direction source P <= 0.0025', limitation='No same-causal-signal, causal-gene or mechanistic replication claim'),
        dict(claim='External validation of shared-concordant architecture', action='HOLD_SIGNAL_LEVEL_CLAIM', evidence='Two Graves representatives supported; no AIHT primary association support', limitation='AIHT nonsignificance is not absence or failed replication'),
        dict(claim='External validation of shared-opposite architecture', action='HOLD_SIGNAL_LEVEL_CLAIM', evidence='GSU056 nonsignificant; GSU254 unavailable in primary set', limitation='No external cross-disease contrast; different signals and ancestry LD not excluded'),
        dict(claim='TSHR Graves predominance / disease specificity', action='DOWNGRADE_EXTERNAL_CLAIM_TO_GRAVES_VARIANT_ASSOCIATION', evidence='GSU231 Graves association supported; AIHT small-effect equivalence not established', limitation='Do not change the original discovery classification'),
        dict(claim='NRG1 Graves predominance', action='HOLD_EXTERNAL_VALIDATION', evidence='GSU138 does not reach either phenotype association threshold', limitation='Nonsignificance does not refute the discovery classification'),
        dict(claim='Absent AIHT effects or definitive nonreplication', action='DO_NOT_CLAIM', evidence='Five evaluable AIHT targets low-precision under OR1.10 scenario; five missing', limitation='No equivalence pass; no absence proof'),
        dict(claim='Cross-ancestry effect equality', action='DO_NOT_CLAIM', evidence='No multiplicity-adjusted heterogeneity positive', limitation='Nonsignificant Q does not establish equality'),
        dict(claim='Two independent BBJ cohort replications or joint disease contrast', action='DO_NOT_CLAIM', evidence='BBJ disease sources share a cohort; covariance unknown', limitation='No pooled meta-analysis or cross-disease difference test'),
        dict(claim='Confirmed divergent-fate mechanism or treatment target', action='DO_NOT_UPGRADE', evidence='Marginal variant associations only; earlier T12/T13/T14 HOLDs unchanged', limitation='No causal or therapeutic validation; manuscript not edited')]
    table(out / 'claim_disposition.tsv', claims)
    forest = []
    z = NormalDist().inv_cdf(.975)
    for r in rows:
        for cohort, prefix in [('FinnGen', 'discovery'), ('BBJ', 'external')]:
            b, s = r[prefix + '_beta'], r[prefix + '_se']
            present = b != 'NA' and s != 'NA'
            if present:
                b, s = float(b), float(s)
                assert math.isfinite(b) and math.isfinite(s) and s > 0
            forest.append(dict(signal_unit_id=r['signal_unit_id'], trait=r['trait'], cohort=cohort,
                discovery_class=r['five_class_discovery'], variant_id=r['variant_id'], nearest_genes=r['nearest_genes'],
                effect_allele=r['effect_allele'], other_allele=r['other_allele'], available=present,
                beta=b, se=s, odds_ratio=math.exp(b) if present else 'NA',
                lower95=math.exp(b - z * s) if present else 'NA', upper95=math.exp(b + z * s) if present else 'NA',
                source_p=r['discovery_p'] if cohort == 'FinnGen' else r['external_source_p'], external_classification=r['classification']))
    table(out / 'forest_data.tsv', forest)
    with (out / 'forest_render.log').open('x') as log:
        subprocess.run(['Rscript', str(ROOT / 'codes/64_plot_t18_forest.R'), str(out)], stdout=log, stderr=subprocess.STDOUT, check=True)
    assert before == {str(p.relative_to(ROOT)): sha(p) for p in inputs}, 'Inputs changed during synthesis'
    assert sha(out / 'priority_variant_validation.tsv') == before[str((SOURCE / 'priority_variant_validation.tsv').relative_to(ROOT))]
    testable = sum(r['availability'] == 'aligned_exact_target' for r in rows)
    report = ['# T18 有限变异层面综合', '', '日期：2026-09-04。完成当前可行的类别汇总、反向候选审查、森林图和证据表述台账；完整 T18 信号级总判定仍 HOLD。', '',
        '## 类别汇总', '', '下表的类别来自 FinnGen 发现集，列出的阳性仅为同病种、同向的代表变异关联支持，不是该类别的外部复制率。分母保持原计划 20 项，不因缺失或分组而重做校正。', '',
        '| 发现类别 | 病种 | 预定 | 可检验 | 同向关联支持 | 未定／低精度 | 不可检验 |', '|---|---|---:|---:|---:|---:|---:|']
    for s in summary:
        if s['planned_variant_tests']:
            report.append('| {discovery_class} | {trait} | {planned_variant_tests} | {testable_variant_tests} | {same_direction_association_support} | {inconclusive_low_precision} | {untestable} |'.format(**s))
    report += ['', 'AIHT-predominant 在发现集中为 0；294 个 unresolved 未纳入本轮，不把未分析记作未复制。GSU058 是原 11 行 priority 文件中的 unresolved 行，不在冻结的 10 个分类明确主候选中。', '',
        '## 可保留的发现与不可升级的结论', '',
        '- Graves 的 GSU052、GSU073（最近基因 LPP）、GSU231（最近基因 TSHR）有同向变异关联支持；最近基因不是因果基因指认。',
        '- T07 原发现集分类保留，但 shared-concordant、shared-opposite、Graves-predominant 均未获外部信号级重验。TSHR 的外部证据只支持 Graves 边际关联，不能证明疾病偏向。',
        '- AIHT 的 5 项可检验结果均未达关联阈值，且固定 OR1.10 情景精度不足；另 5 项不可检验。不能说无效应或未复制。7 项低精度分类不代表只有 7 项情景功效不足，3 项关联阳性在该固定情景下也低于 80%。',
        '- 各按 20 项校正后 Q 阳性为 0、TOST 等效通过为 0；前者不证明效应一致，后者不证明另一病近零，也不证明效应超出等效界限。', '',
        '## shared-opposite 严格审查', '',
        '- GSU056：两病代表变异均可检验，源行、双向坐标和等位基因方向已有全量核验；BBJ 两病估计均为负，仅为描述性符号。Graves P≈0.0513、AIHT P≈0.6718，均未达 0.0025。AIHT 与发现集相反的符号不是反向关联支持；其 Q 未校正 P≈0.049、校正 P≈0.980，不是显著异质性。',
        '- GSU254：回文代表变异不在冻结主分析集，两病均不可检验；不以代理 SNP 补位。',
        '- 协调复核不等于排除不同信号或祖源 LD 差异；缺少合格日本／东亚 LD，无法完成最严格的信号级排除。BBJ 两病种协方差未知，不检验两病效应差异；现有结果既不验证也不推翻 shared-opposite 类别。', '',
        '## 森林图与统计边界', '',
        '每个病种单列一幅图：10 个预定代表位点均保留，FinnGen 与 BBJ 按相同效应等位基因展示 OR 和普通正态近似 95% CI；BBJ 缺失不补零。两幅使用相同横轴范围，不跨病种作效应差异推断、不合并效应。',
        '关联支持沿用源双侧 P≤0.0025，不依据图中 95% CI 是否跨 1 判定。普通 Wald CI、Q、TOST、精度情景不继承 SAIGE/SPA 校准；两队列表型并非已证实等同。方法来源沿用 [冻结验证报告](../../../documentation/analysis/T16_T17_FEASIBLE_VARIANT_VALIDATION.md)。', '',
        '![Graves variant effects](forest_Graves.png)', '', '![AIHT variant effects](forest_AIHT.png)', '',
        '## 验证、文件与停止边界', '',
        '输入清单记录当前已验证结果边界、冻结合同和脚本 SHA-256，执行前后保持一致；本轮不重复扫描多 GB 原始数据或协调数据库。复制统计表逐字节不变，在新目录用独立 base R 再算 29/29 项通过。',
        '`class_summary.tsv` 为类别计数；`representative_evidence.tsv` 保留逐位点双病种证据；`shared_opposite_audit.tsv` 保留未排除事项；`claim_disposition.tsv` 是分析表述台账，不代表已修改稿件。`forest_data.tsv` 是可编辑图数据，SVG 为矢量导出，PNG 为审阅版。',
        '尚需：祖源匹配 LD 和信号级分析、表型可比性及独立欧洲来源。完整 T18 保持 HOLD，T12/T13/T14 原 HOLD 不变。本轮未改变 T07 分类、阈值、原统计结果或稿件，不进入 T19。', '']
    (out / 'T18_LIMITED_SYNTHESIS.md').write_text('\n'.join(report))
    dump(out / 'decision.json', dict(status='PASS_T18_VARIANT_SYNTHESIS_HOLD_SIGNAL_LEVEL',
        completed_at_utc=datetime.now(timezone.utc).isoformat(), planned_tests=20, testable_tests=testable,
        variant_level_synthesis_completed=True, independent_numeric_checks_passed=len(checks),
        input_hashes_unchanged=True, forest_plots_generated=True, visual_audit='PENDING_INDEPENDENT_REVIEW',
        signal_level_validation_completed=False, shared_opposite_strict_signal_audit_completed=False,
        cross_disease_difference_executed=False, manuscript_changed=False, t19_started=False,
        full_t18_complete=False))
    print(json.dumps({'output': str(out), 'testable': testable, 'R_checks': len(checks), 'status': 'NUMERICAL_PASS_VISUAL_REVIEW_PENDING'}))


if __name__ == '__main__':
    main()
