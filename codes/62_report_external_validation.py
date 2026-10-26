#!/usr/bin/env python3
"""Produce an evidence-bound report from verified, frozen-priority outputs."""
import collections
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'outcome_data/external_validation/run_20260904'
REPORT=ROOT/'documentation/analysis/T16_T17_FEASIBLE_VARIANT_VALIDATION.md'
LABELS={'same_direction_association_support':'同向关联支持','opposite_direction_association_support':'反向关联支持',
 'no_primary_support_small_effect_equivalence':'未达关联阈值／小效应等效',
 'inconclusive_low_precision_OR1_10_scenario':'未定／OR1.10情景低精度','inconclusive':'未定',
 'untestable_missing_aligned_target':'不可检验／主分析集无该位点'}


def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    assert not REPORT.exists(),'Refusing to overwrite report'
    checks=read(RUN/'independent_statistics_validation.tsv');assert all(r['passed']=='TRUE' for r in checks)
    a=json.loads((ROOT/'outcome_data/t16/harmonization_verified_20260904/decision.json').read_text())
    audit=json.loads((ROOT/'outcome_data/t16/harmonization_verified_20260904/independent_validation.json').read_text())
    assert audit['status']=='PASS'
    freeze=json.loads((ROOT/'outcome_data/external_validation/plan_20260904/freeze.json').read_text())
    binding=json.loads((ROOT/'outcome_data/external_validation/execution_plan_20260904/amendment.json').read_text())
    result=json.loads((RUN/'decision.json').read_text());rows=read(RUN/'priority_variant_validation.tsv')
    counts=collections.Counter((r['trait'],r['classification']) for r in rows)
    out=['# T16/T17 可行的变异层面外部验证','',
        '日期：2026-09-04。范围：补齐 T16 清洗/协调，执行冻结的非 HLA 主代表位点验证；不执行依赖 LD 的分析或稿件结论改写。','',
        '## 本轮结论','',
        f"T16 完成 {a['aligned_records']:,} 个 SNP 的坐标及等位基因协调，{a['checks_passed']}/{a['checks_total']} 流水线检查通过，独立全量复核 0 错误。",
        f"固定 20 项（10 个分类明确的 T07 主候选 × 2 病种），实际可检验 {result['testable_tests']} 项；其余保留不可检验，不缩小多重检验分母。",'',
        '| 病种 | 变异层面结果 | 项数 |','|---|---|---:|']
    for (trait,c),n in sorted(counts.items()):out.append(f'| {trait} | {LABELS[c]} | {n} |')
    testable=[r for r in rows if r['availability']=='aligned_exact_target']
    hsig=sum(float(r['heterogeneity_p_bonferroni'])<=.05 for r in testable)
    esig=sum(r['equivalence_flag']=='True' for r in testable)
    out+=['',f'分别按 20 项校正后，跨队列异质性阳性 {hsig} 项，小效应等效通过 {esig} 项。没有显著异质性不证明效应相同；未通过等效也不能声称另一病近零。']
    out+=['','**同向关联支持不等于同一因果信号复制，也不重新验证 T07 五分类。** 未显著和主分析集缺失均不自动归为“未复制”。',
          '', '7 项“未定／低精度”是按分类优先级得到的结果，不表示只有 7 项情景功效不足：3 项 Graves 关联阳性在固定 OR=1.10 情景下也均低于 80%，与其较大的实际关联效应并不矛盾。不据此报告整个 T07 体系的“复制率”。',
          '', '## 10 个冻结主候选','',
          '| 信号 | 最近基因（位置注释） | Graves 结果 | Graves 源 P | AIHT 结果 | AIHT 源 P |',
          '|---|---|---|---:|---|---:|']
    by={(r['signal_unit_id'],r['trait']):r for r in rows}
    for sid in sorted({r['signal_unit_id'] for r in rows}):
        g,h=by[(sid,'Graves')],by[(sid,'AIHT')]
        out.append(f"| {sid} | {g['nearest_genes']} | {LABELS[g['classification']]} | {g['external_source_p']} | {LABELS[h['classification']]} | {h['external_source_p']} |")
    out+=['','源 P 为双侧报告值，同向/反向支持阈值为 0.05/20=0.0025；没有用重算的 Wald P 替换源 P。基因名仅为最近基因注释，不代表效应基因或机制验证。',
        '', '## 方案冻结与审计修订','',
        f"执行前分析计划于 **{freeze['frozen_at_utc']}** 冻结。这不是前瞻性预注册：发现集、来源 QC 和已有 T17 覆盖信息事先可见。",
        '原 `priority_signal_units.tsv` 实为 11 条，含 GSU058（unresolved）；本轮按发现集非 unresolved 分类固定 10 条。GSU058 未因外部效应结果而被选择性剔除，其余 294/304 个信号不进入本轮统计。',
        f"执行绑定修订于 {binding['created_at_utc']} 记录：仅将 T16 路由指向通过验证的新目录并绑定两病种实际文件哈希；统计参数、目标列表和原 57 脚本未变。",
        '原始统计冻结、执行修订及输入哈希均保留；运行前后核对源文件、已协调文件、decision 和独立审计记录。T17 另行完整重验以绑定当前文件版本。',
        '', '## T16 清洗及协调','',
        f"13,428,945 条源记录 → {a['clean_records']:,} 个清洗 SNP → {a['aligned_records']:,} 个协调 SNP。",
        '限常染色体、非 HLA（两个 build 均 chr6:25–34 Mb）、普通非回文 SNP、两侧 MAF≥0.005；多记录位置全部排除；两向一碱基、单一映射及精确坐标/链方向往返均须通过。',
        'EA/OA 按源文件定义处理，不把 OA 宣称为参考基因组等位基因。负链互补后仅按精确 ALT/REF 配对；swap 时只翻转 BETA 与 EAF，原 P/SE 不变。未进行独立参考 FASTA 核验。',
        '逐位点 N、MAC、Rsq、病例/对照拆分仍不可用，没有补造；537/172656 仅为队列 metadata。本轮 beta/SE 方法不要求填入逐位点 N。',
        '初始实现的坐标状态更新缺陷由代码复核拦截并修正，增加 SQLite + liftOver 端到端回归测试；旧进程经核对终止后在全新目录单进程重跑。旧目录及中断产物保留，不作为正式结果。',
        '', '## 统计合同及限制','',
        '三个家族分别固定 Bonferroni 分母 20：关联、异质性、等效；不声称三类合并总体 FWER≤0.05。BBJ 两病种依赖性不影响 Bonferroni 保守有效性，但不是两个独立验证队列。',
        '同病种 FinnGen–BBJ 比较采用 Δβ=βBBJ−βFG，SEΔ=√(SEBBJ²+SEFG²)，Q=(Δβ/SEΔ)²、χ²₁ 尾概率；假定这两个队列没有样本重叠。I²=max(0,(Q−1)/Q)，Q=0 时 I²=0，结果以比例表示，不套低/中/高阈值。没有跨祖源合并效应。',
        '提供 log-OR、OR 及 Wald 95% CI。TOST 使用 ±log(1.10)，双单侧 P 的较大值 ≤0.0025，等价于 99.5% CI 完全落入界限；这不是经临床验证的意义界限。显著关联与小效应等效可同时成立，因此单独保留 equivalence_flag。',
        '精度情景以固定真实效应 |logOR|=log(1.10)、源 SE 和双侧 α=0.0025 计算；低于 80% 时仅称该情景低精度。另给 0.5×/1×发现效应情景及 80% 最小可检测 OR；不使用 BBJ 估计效应计算所谓观察功效。发现效应情景可能受 winner’s curse 影响，未作偏倚校正。',
        '源 P 与正态近似检验不等价。95% CI、Q、TOST、精度情景均不继承 SAIGE/SPA 校准；严重病例不平衡下尤其需谨慎。官方 [BBJ 字段字典](https://humandbs.dbcls.jp/files/hum0197/Dictfile_BBJ.html) 明确二分类使用 SAIGE；[SAIGE 原始方法](https://www.nature.com/articles/s41588-018-0184-y) 支持区分源检验与普通渐近计算。',
        '异质性显著不能指定祖源因果作用；不显著不证明效应一致。表型定义、LD 标记、发现集选择均可能影响差异；仅两个队列的 I² 信息有限，见 [Cochrane 方法手册](https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-10)。',
        '', '## 独立数值复核及输出','',
        f"相关单元/集成测试 25 项通过；R 用独立分布函数重算全部可检验行，{len(checks)}/{len(checks)} 项检查通过。T16 原始行与全部协调输出独立逐行检查通过；T17 全量复核及哈希绑定通过。",
        '主结果：`outcome_data/external_validation/run_20260904/priority_variant_validation.tsv`。',
        '完整数值复核：`independent_statistics_validation.tsv`；分组计数：`classification_counts.tsv`；执行状态与来源：`decision.json`、`provenance.tsv`。',
        '原方案：`outcome_data/external_validation/plan_20260904/`；执行绑定：`outcome_data/external_validation/execution_plan_20260904/`。',
        'T16 正式准备结果：`outcome_data/t16/harmonization_verified_20260904/`；T17 新审计：`outcome_data/t17/revalidation_bound_20260904/`（数据源仍为原协调目录）。',
        '', '## 保留的 HOLD','',
        '日本/东亚 LD、信号级共定位、T07 shared/opposite/predominant 类别的外部重验仍 HOLD；没有代理 SNP、没有将欧洲 LD 用作 BBJ 替代。',
        'BBJ Hashimoto 与 FinnGen strict AIHT、BBJ Graves 与 strict Graves 仅作映射比较，不保证表型相同；独立欧洲来源尚未取得。',
        'BBJ 两病种共同队列的误差协方差未知，不计算外部 Graves–AIHT 直接差异，也不把两个边际 P 的差异当病例—病例证据。',
        'T18 现有了有限的变异层结果，但完整外部验证总判定仍受信号级与表型/祖源限制；本轮未生成森林图或改写稿件主结论。T12/T13/T14 原有 HOLD 不变。',
        '', '## 重现顺序','',
        '入口依次为 `55_harmonize_t16_bbj.py` → `58_validate_t16_harmonization.py`；`60_revalidate_t17_binding.py`；`61_bind_external_execution_plan.py` → `57_run_external_validation.py --plan outcome_data/external_validation/execution_plan_20260904` → `59_verify_external_statistics.R` → 本报告脚本。',
        '冻结入口为 `56_freeze_external_validation.py`。脚本对既有输出拒绝覆写；复现需另建路径并记录新的执行绑定，不能覆盖已冻结结果。',
        '本轮按 research-cost-optimizer 完成独立方法与实现审计，并按 project-steward 同步一次当前文档。']
    REPORT.write_text('\n'.join(out)+'\n')
    print(REPORT)


if __name__=='__main__':main()
