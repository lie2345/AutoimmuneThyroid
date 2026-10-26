#!/usr/bin/env python3
"""Build publication views of frozen results; never rerun statistical analyses."""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

REV = Path(__file__).resolve().parents[1]
ROOT = REV.parents[3]
BASE = ROOT / "outputs/sci_writing_20261017/revisions/circular_v1"
T19 = ROOT / "outputs/t19_t24_20261017/t19/tables"
T20 = ROOT / "outputs/t19_t24_20261017/t20"
TABLES, DATA, REVIEW = (REV / p for p in ("tables", "source_data", "reviews"))
for directory in (TABLES, DATA, REVIEW):
    directory.mkdir(parents=True, exist_ok=True)


def rel(path):
    return str(Path(path).relative_to(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


@dataclass
class Source:
    path: Path
    line: int
    data: dict


inputs = {}


def read(path):
    path = Path(path)
    inputs[rel(path)] = sha(path)
    records = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        _ = reader.fieldnames
        last = reader.line_num
        for row in reader:
            records.append(Source(path, last + 1, row))
            last = reader.line_num
    return records


def tsv(path, headers, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


provenance = []


def cell(table, row, column, template, *refs):
    values = [source.data[key] for source, key in refs]
    value = template.format(*values)
    for (source, key), original in zip(refs, values):
        provenance.append({"table": table, "row_id": row, "column": column,
                           "source_file": rel(source.path), "source_line": source.line,
                           "source_column": key, "source_value": original,
                           "display_value": value, "relation": "verbatim field in publication cell"})
        assert original in value, (table, row, column, key)
    return value


def payload(name, ident, title, columns, rows, notes):
    result = {"id": ident, "title": title,
              "columns": [{"key": key, "label": label} for key, label in columns],
              "rows": rows, "notes": notes, "references": []}
    dump(TABLES / f"{name}.json", result)
    tsv(TABLES / f"{name}.tsv", ["row_id"] + [key for key, _ in columns], rows)
    return result


def retain_full(name, records):
    source = records[0].path
    target = DATA / f"{name}.tsv"
    shutil.copyfile(source, target)
    assert source.read_bytes() == target.read_bytes()
    dump(DATA / f"{name}.json", {"source_file": rel(source), "source_sha256": sha(source),
                               "columns": list(records[0].data), "rows": [r.data for r in records]})


summary = read(T19 / "hla_within_trait_summary.tsv")
flow = read(T19 / "hla_qc_filter_flow.tsv")
leads = read(T19 / "hla_regional_lead_variants.tsv")
cross = read(T19 / "hla_cross_trait_at_regional_leads.tsv")
endpoints = read(T20 / "endpoint_metadata.tsv")
panel = read(T20 / "frozen_panel.tsv")
tests = read(T20 / "frozen_tests.tsv")
lookup = read(T20 / "lookup_results.tsv")
freeze = read(T20 / "freeze_record.tsv")
confirm = read(T20 / "hla_selection_confirmation.tsv")
priority = read(ROOT / "outcome_data/external_validation/run_20260904/priority_variant_validation.tsv")

for name, records in [("table_s7a_region_full", summary), ("table_s7a_qc_full", flow),
                      ("table_s7b_full", leads), ("table_s7c_full", cross),
                      ("table_s8a_full", endpoints), ("table_s8_frozen_panel", panel),
                      ("table_s8_frozen_tests", tests), ("table_s8_lookup_results", lookup),
                      ("table_s8_freeze_record", freeze), ("table_s8_hla_selection_confirmation", confirm)]:
    retain_full(name, records)

# Main-table row content is inherited; only the two approved endpoint rows are added.
main = []
for name in ("table1", "table2"):
    source = BASE / "tables" / f"{name}.json"
    inputs[rel(source)] = sha(source)
    main.append(json.loads(source.read_text(encoding="utf-8")))
table1, table2 = main
original1 = list(table1["rows"])
original2 = json.loads(json.dumps(table2["rows"]))
for source in endpoints:
    d = source.data
    endpoint = d["endpoint"]
    rid = "FG_GO" if endpoint == "GRAVES_OPHT" else "FG_GO_STRICT"
    label = "Graves ophthalmopathy" if rid == "FG_GO" else "Graves ophthalmopathy, strict"
    row = {"row_id": rid, "section": "Supporting evidence",
           "data": cell("table1", rid, "data", "FinnGen R13 " + label + "\n{0}", (source, "endpoint")),
           "definition": cell("table1", rid, "definition", "{0}. Full endpoint rules are given in Supplementary Table S8A.", (source, "case_algorithm_source_summary")),
           "sample_unit": cell("table1", rid, "sample_unit", "{0} cases; {1} controls", (source, "n_cases"), (source, "n_controls")),
           "role_limit": "FinnGen R13; GRCh38; population noncase controls. Exploratory lookup of variants selected after discovery. Individual overlap is unquantified; this is not independent replication."}
    table1["rows"].append(row)
table1["notes"].extend([
    "The descriptive HLA analysis reuses the two primary FinnGen R13 GWAS datasets over GRCh38 chr6:25,000,000–34,000,000 inclusive (Supplementary Table S7). It does not add a cohort or estimate independent HLA signals.",
    "The ophthalmopathy endpoints use population noncase controls: ‘Controls for this endpoint are individuals that are not cases.’ They do not compare ophthalmopathy cases with Graves patients without ophthalmopathy. Individual overlap and strict nesting of the two endpoint case sets cannot be established from these summary data. The Graves discovery endpoint itself includes GRAVES_OPHT in its INCLUDE rule; the exploratory lookup is not independent replication (Supplementary Table S8)."
])
table2["notes"].append("Exploratory ophthalmopathy lookups of the existing ten non-HLA representatives and two separately selected HLA regional leads are reported in Supplementary Table S8. These population-level associations provide phenotype context and do not validate component identity, H4/H3 classifications or ophthalmopathy modification conditional on Graves’ disease.")
assert table1["rows"][:len(original1)] == original1
assert table2["rows"] == original2 and len(original2) == 10
for name, data in zip(("table1", "table2"), main):
    dump(TABLES / f"{name}.json", data)
    columns = ["row_id"] + [column["key"] for column in data["columns"]]
    tsv(TABLES / f"{name}.tsv", columns, data["rows"])
    shutil.copyfile(BASE / "source_data" / f"{name}_source_map.tsv", DATA / f"{name}_source_map.tsv")
shutil.copyfile(BASE / "tables/references.json", TABLES / "references.json")

hla_note = "HLA is GRCh38 chr6:25,000,000–34,000,000 inclusive. Counts describe retained variant rows, not independent signals. The primary GWAS samples are reused; no HLA cohort was added."
ci_note = "OR and approximate normal 95% confidence limits are copied exactly from the frozen outputs (exp[beta ± 1.959963984540054 × SE]). Original source P values are retained; these intervals are not inversions of the source GWAS calibrated test. No model was fitted or P value reconstructed for this table."
lead_note = "Regional leads minimize positive source P among variants retained by the established general QC. All exact minimum-P ties are retained. Different leads do not establish distinct causal signals; same-direction effects do not establish a shared causal signal. No conditional HLA analysis or causal-gene assignment was performed."

region_rows = []
for s in summary:
    rid = s.data["trait_id"]
    C = lambda col, template, *keys: cell("table_s7a_region", rid, col, template, *[(s,k) for k in keys])
    region_rows.append({"row_id": rid,
        "trait": C("trait", "{0}\n{1}", "trait_id", "label"),
        "region": C("region", "{0}; chr{1}:{2}–{3} inclusive", "genome_build", "hla_chr", "hla_start_inclusive", "hla_end_inclusive"),
        "samples": C("samples", "Cases {0}\nControls {1}", "n_cases", "n_controls"),
        "retained": C("retained", "Raw {0}\nRetained {1}\nMAF ≥ {2}", "raw_hla_rows", "retained_hla_rows", "general_maf_min"),
        "variant_types": C("variant_types", "SNP {0}\nIndel/other {1}\nPalindromic {2}", "retained_snp_rows", "retained_indel_or_other_rows", "retained_palindromic_rows"),
        "source_p_counts": C("source_p_counts", "P < {0}: {1}\nP < {2}: {3}", "genomewide_threshold", "source_p_lt_genomewide", "suggestive_threshold", "source_p_lt_suggestive"),
        "extent": C("extent", "Retained positions {0}–{1}\nRows at start/end {2}/{3}\nRaw P=0 rows {4}", "retained_min_pos", "retained_max_pos", "rows_at_start_boundary", "rows_at_end_boundary", "raw_source_p_zero_rows")})
supp = []
supp.append(payload("table_s7a_region", "Table S7A", "HLA region and retained variant counts",
    [("trait","Trait"),("region","Region"),("samples","Released GWAS samples"),("retained","Raw and retained rows"),("variant_types","Retained variant types"),("source_p_counts","Source P thresholds"),("extent","Coverage and boundary QC")], region_rows, [hla_note]))
flow_rows=[]
for i,s in enumerate(flow,1):
    rid=f"QC{i:02d}"
    row={"row_id":rid}
    for key in s.data:
        row[key]=cell("table_s7a_qc",rid,key,"{0}",(s,key))
    flow_rows.append(row)
supp.append(payload("table_s7a_qc", "Table S7A QC", "Sequential HLA quality-control flow",
    [("trait_id","Trait"),("step","QC step"),("rows_before","Rows before"),("rows_after","Rows after"),("removed","Rows removed")], flow_rows,
    [hla_note,"The ordered steps and all counts are copied from the existing QC flow. The MAF threshold is 0.005; duplicate variant IDs use minimum source P. Detailed inclusion rules are retained in the frozen source table."]))

lead_rows=[]
for s in leads:
    rid=s.data["trait_id"]
    C=lambda col,template,*keys:cell("table_s7b",rid,col,template,*[(s,k) for k in keys])
    lead_rows.append({"row_id":rid,
        "selected_for_trait":C("selected_for_trait","{0}","trait_id"),
        "variant":C("variant","{0}\n{1}\nREF/ALT {2}/{3}\nEffect/other {4}/{5}","variant_id","rsid","ref","alt","effect_allele","other_allele"),
        "frequency":C("frequency","EAF {0}\nMAF {1}","eaf","maf"),
        "beta_se":C("beta_se","Beta {0}\nSE {1}","beta","se"),
        "or_ci":C("or_ci","OR {0}\n95% CI {1}–{2}","odds_ratio","or_ci95_lower","or_ci95_upper"),
        "source_p":C("source_p","P {0}\nSource −log10(P) {1}\n{2}","source_p","source_mlogp","source_p_status"),
        "selection":C("selection","{0}\nExact minimum-P ties {1}\nSNP {2}; palindromic {3}\nCases/controls {4}/{5}","lead_selection","min_p_ties","is_snp","palindromic","n_cases","n_controls")})
supp.append(payload("table_s7b","Table S7B","Trait-selected HLA regional leads",
    [("selected_for_trait","Selected for trait"),("variant","Exact variant and alleles"),("frequency","Allele frequencies"),("beta_se","Marginal effect"),("or_ci","OR and approximate normal 95% CI"),("source_p","Source P and source mlogp"),("selection","Selection and samples")],lead_rows,[hla_note,lead_note,ci_note]))

cross_rows=[]
for s in cross:
    rid=s.data["selected_as_lead_for"]+"__"+s.data["evaluated_trait"]
    C=lambda col,template,*keys:cell("table_s7c",rid,col,template,*[(s,k) for k in keys])
    cross_rows.append({"row_id":rid,
        "selected_for_trait":C("selected_for_trait","{0}","selected_as_lead_for"),
        "evaluated_trait":C("evaluated_trait","{0}\nCases/controls {1}/{2}","evaluated_trait","n_cases","n_controls"),
        "variant":C("variant","{0}\n{1}\nEffect/other {2}/{3}","variant_id","rsid","effect_allele","other_allele"),
        "frequency":C("frequency","EAF {0}\nMAF {1}","eaf","maf"),
        "beta_se":C("beta_se","Beta {0}\nSE {1}","beta","se"),
        "or_ci":C("or_ci","OR {0}\n95% CI {1}–{2}","odds_ratio","or_ci95_lower","or_ci95_upper"),
        "source_p_matching":C("source_p_matching","P {0}\nSource −log10(P) {1}\n{2}\nDirection: {3}","source_p","source_mlogp","matching_status","effect_direction_relative_to_selecting_trait")})
supp.append(payload("table_s7c","Table S7C","Same-variant marginal effects at both HLA leads",
    [("selected_for_trait","Lead selected for"),("evaluated_trait","Evaluated trait and samples"),("variant","Variant and effect allele"),("frequency","Allele frequencies"),("beta_se","Marginal effect"),("or_ci","OR and approximate normal 95% CI"),("source_p_matching","Source P and matching")],cross_rows,[hla_note,lead_note,ci_note,"All four same-variant estimates are shown. No cross-trait effect-difference test was performed."]))

endpoint_notes=["Both endpoints use population noncase controls, rather than Graves patients without ophthalmopathy. Exact individual overlap with discovery or with the other ophthalmopathy endpoint is unavailable. Strict nesting of the case sets is unconfirmed.","Variant selection occurred after the original Graves–AIHT discovery. All endpoints are FinnGen R13; the lookup is not independent replication. The Graves discovery endpoint includes GRAVES_OPHT in its INCLUDE rule."]
endpoint_rows=[]
for s in endpoints:
    rid=s.data["endpoint"]
    C=lambda col,template,*keys:cell("table_s8a",rid,col,template,*[(s,k) for k in keys])
    name="Graves ophthalmopathy" if rid=="GRAVES_OPHT" else "Graves ophthalmopathy, strict"
    endpoint_rows.append({"row_id":rid,
        "endpoint":C("endpoint",name+"\n{0}","endpoint"),
        "release_build":C("release_build","{0}\n{1}","release","genome_build"),
        "counts":C("counts","Cases {0}\nPopulation noncase controls {1}","n_cases","n_controls"),
        "case_rule":C("case_rule","{0}","case_algorithm_source_summary"),
        "control_rule":C("control_rule","{0}\n{1}","control_definition_verbatim","control_definition_status"),
        "definition_source":C("definition_source","{0}\nSnapshot {1}, control-rule line {2}\nCount source {3}, lines {4}","official_url","definition_snapshot","definition_line","count_source","count_line")})
supp.append(payload("table_s8a","Table S8A","FinnGen R13 ophthalmopathy endpoint definitions",
    [("endpoint","Official endpoint"),("release_build","Release and build"),("counts","Released GWAS samples"),("case_rule","Case algorithm from official snapshot"),("control_rule","Verbatim control definition"),("definition_source","Definition and count sources")],endpoint_rows,endpoint_notes))

panel_by_id={s.data["variant_id"]:s for s in panel}
test_by_key={(s.data["endpoint"],s.data["variant_id"]):s for s in tests}
gsu={}
for s in priority:
    variant=s.data["variant_id"]
    previous=gsu.setdefault(variant,s.data["signal_unit_id"])
    assert previous==s.data["signal_unit_id"]
freeze_data=freeze[0].data
assert freeze_data["panel_sha256"]==sha(T20/"frozen_panel.tsv")
assert freeze_data["tests_sha256"]==sha(T20/"frozen_tests.tsv")
assert freeze_data["source_priority_sha256"]==sha(priority[0].path)
assert len(lookup)==24 and len(panel)==12 and len(leads)==2 and len(cross)==4
assert len(test_by_key)==24

all_full=[]
compact={"non_hla":[],"hla":[]}
for s in lookup:
    d=s.data
    p=panel_by_id[d["variant_id"]]
    t=test_by_key[(d["endpoint"],d["variant_id"])]
    for key,value in p.data.items():
        assert d[key]==value,(d["variant_id"],key)
    group="non_hla" if int(d["panel_index"])<=10 else "hla"
    selected_id=gsu[d["variant_id"]] if group=="non_hla" else d["selection_class"]
    full=dict(d)
    full.update({"panel_group":"Non-HLA discovery representatives" if group=="non_hla" else "HLA regional leads",
                 "historical_signal_unit_id":gsu[d["variant_id"]] if group=="non_hla" else "NA",
                 "test_index":t.data["test_index"],"bonferroni_family_n":t.data["bonferroni_family_n"],
                 "lookup_source_file":rel(s.path),"lookup_source_line":str(s.line),
                 "panel_source_file":rel(p.path),"panel_source_line":str(p.line),
                 "test_source_file":rel(t.path),"test_source_line":str(t.line),
                 "frozen_at_utc":freeze_data["frozen_at_utc"],"frozen_panel_sha256":freeze_data["panel_sha256"],
                 "frozen_tests_sha256":freeze_data["tests_sha256"],
                 "selection_confirmation":"Historical frozen representative" if group=="non_hla" else "T19 lead identity confirmed after freeze; frozen panel unmodified"})
    all_full.append(full)
    rid="T"+t.data["test_index"].zfill(2)
    table="table_s8b_"+group
    C=lambda col,template,*keys:cell(table,rid,col,template,*[(s,k) for k in keys])
    compact[group].append({"row_id":rid,
        "selection_variant":C("selection_variant","Panel {0}; "+selected_id+"\n{1}\n{2}","panel_index","variant_id","rsids_source"),
        "endpoint_samples":C("endpoint_samples","{0}\nCases/controls {1}/{2}","endpoint","n_cases","n_controls"),
        "alleles_frequency":C("alleles_frequency","REF/ALT {0}/{1}\nEffect/other {2}/{3}\nEAF_ALT {4}; MAF {5}","ref","alt","effect_allele","other_allele","eaf_source_string","maf"),
        "match_qc":C("match_qc","{0}\n{1}; {2}\nPosition/exact rows {3}/{4}","match_status","qc_status","qc_reason","source_position_row_count","exact_ref_alt_row_count"),
        "beta_se":C("beta_se","Beta {0}\nSE {1}","beta_source_string","se_source_string"),
        "or_ci":C("or_ci","OR {0}\n95% CI {1}–{2}","OR_ALT","OR_lower95_normal_approx","OR_upper95_normal_approx"),
        "source_p":C("source_p","P {0}\nOptional min(1,24×P) {1}\nP < 0.05/24: {2}","source_p_string","p_bonferroni24_optional_descriptive","source_p_below_0_05_over_24"),
        "direction_context":C("direction_context","Discovery class: {0}\nDirection vs Graves: {1}\nDirection vs AIHT: {2}","five_class_discovery","direction_vs_graves_discovery","direction_vs_aiht_discovery")})
full_columns=list(all_full[0])
tsv(DATA/"table_s8b_all24_full.tsv",full_columns,all_full)
dump(DATA/"table_s8b_all24_full.json",{"id":"Table S8B full machine source","columns":full_columns,"rows":all_full})
freeze_note=("The panel was frozen at "+freeze_data["frozen_at_utc"]+" before this run’s eye-data extraction. It contains 10 historical non-HLA discovery representatives and two HLA marginal regional leads, crossed with two endpoints: the complete family is 24 lookups. Panel SHA-256: "+freeze_data["panel_sha256"]+". Frozen tests SHA-256: "+freeze_data["tests_sha256"]+".")
lookup_notes=endpoint_notes+[ci_note,freeze_note,
    "All 24 planned endpoint–variant combinations are retained, including results above nominal or optional thresholds. The optional adjustment is the saved min(1,24×source P), with the full family of 24. It is descriptive after discovery selection, not confirmatory validation. Source P values are not recomputed from beta/SE.",
    "Non-HLA GSU identifiers refer to the historical discovery representatives in the frozen priority source. The HLA selection source was primary-agent assignment at freeze; the two lead identities were subsequently confirmed against T19 without changing the frozen panel. Full selection sources, confirmation timestamps, query, matching, QC, missingness, discovery betas and all original fields accompany the machine-readable tables.",
    "REF and ALT refer to GRCh38. Effects are aligned to ALT; EAF_ALT is its frequency. Discovery direction and class are historical context only. Nearest-gene fields are retained only in the full source tables and do not identify causal genes. NA in historical class/direction for HLA means unavailable or inapplicable; it is not a negative effect."]
lookup_columns=[("selection_variant","Frozen selection and exact variant"),("endpoint_samples","Endpoint and samples"),("alleles_frequency","Alleles and frequencies"),("match_qc","Matching and general QC"),("beta_se","Marginal effect"),("or_ci","OR and approximate normal 95% CI"),("source_p","Source P and optional adjustment"),("direction_context","Historical discovery context")]
for group,label in (("non_hla","Non-HLA discovery representatives"),("hla","HLA regional leads")):
    supp.append(payload("table_s8b_"+group,"Table S8B "+("non-HLA" if group=="non_hla" else "HLA"),label,lookup_columns,compact[group],lookup_notes))

# Compact views carry exact lexical source fields. Full source copies preserve all columns.
all_payloads={name:json.loads((TABLES/f"{name}.json").read_text(encoding="utf-8")) for name in
              ("table1","table2","table_s7a_region","table_s7a_qc","table_s7b","table_s7c","table_s8a","table_s8b_non_hla","table_s8b_hla")}
dump(TABLES/"table_schema_manifest.json",{
    "compact_schema":"{id,title,columns:[{key,label}],rows:[{row_id,<column-key>:value}],notes:[text],references:[object]}",
    "publication_tables":{name:{"rows":len(p["rows"]),"columns":[c["key"] for c in p["columns"]]} for name,p in all_payloads.items()},
    "full_source_tables":"source_data/table_s7*_full.{tsv,json}; source_data/table_s8a_full.{tsv,json}; source_data/table_s8b_all24_full.{tsv,json}",
    "machine_values":"All raw source fields retain source text exactly. Full merged S8B adds provenance, frozen-test fields and historical GSU identifiers without changing any lookup field.",
    "s7_panels":["table_s7a_region","table_s7a_qc","table_s7b","table_s7c"],
    "s8_panels":["table_s8a","table_s8b_non_hla","table_s8b_hla"]})
tsv(DATA/"table_supplementary_source_map.tsv",list(provenance[0]),provenance)
dump(DATA/"table_supplementary_source_map.json",provenance)
# Add endpoint rows to the existing main-table provenance without changing its legacy fields.
map_path=DATA/"table1_source_map.tsv"
legacy=read(map_path)
legacy_fields=list(legacy[0].data)
new_main=[{"table":"table1","row_id":p["row_id"],"column":p["column"],"source_file":p["source_file"],
           "source_locator":"physical line "+str(p["source_line"]),"evidence_key":p["source_column"],"relation":p["relation"]}
          for p in provenance if p["table"]=="table1"]
new_main.extend([
    {"table":"table1","row_id":"note_7","column":"notes","source_file":rel(T19/"hla_within_trait_summary.tsv"),"source_locator":"physical lines 2–3","evidence_key":"reused primary traits and HLA boundaries","relation":"descriptive range; no new cohort"},
    {"table":"table1","row_id":"note_8","column":"notes","source_file":"outputs/manuscript_increment_t19_t20_20261018/sources/strict_graves_risteys_visible.txt","source_locator":"physical lines 79–96","evidence_key":"INCLUDE E4_THYTOXGOITDIF|GRAVES_OPHT","relation":"discovery-phenotype composition"}])
tsv(map_path,legacy_fields,[s.data for s in legacy]+new_main)

notes=["# Publication tables and source boundaries","",
       "The tables present existing frozen outputs. No statistical analysis, association model, fine-mapping or colocalization was rerun.",""]
for data in list(main)+supp:
    notes.extend(["## "+data["id"]+". "+data["title"],""])
    notes.extend(str(i)+". "+note for i,note in enumerate(data["notes"],1))
    notes.append("")
notes.extend(["## Machine-readable detail", "",
    "The compact TSV/JSON views are intended for Word and publication formatting. The source_data directory retains byte-identical full source tables and JSON copies, plus an all-24 merged lookup table containing every original lookup field. The compact values retain the exact source text, including source P strings and CI endpoints. No rounding is introduced.","",
    "Every copied field used in a compact supplementary view has source file, physical line, source column and exact source text in table_supplementary_source_map.tsv. Existing main-table evidence maps and references are retained. Bibliographic claims in Table 2 are inherited from the accepted baseline and were not re-audited against the original publications during this formatting task.","",
    "The workbooks are editable static presentations. Editing workbook cells does not update the frozen analysis. A refresh uses the retained builders and the same authoritative source files."])
(TABLES/"table_notes.md").write_text("\n".join(notes)+"\n",encoding="utf-8")

# Agreed checks: counts, source identity, exact lexical copies, and scope preservation.
for p in provenance:
    source_rows=read(ROOT/p["source_file"])
    original=next(s.data[p["source_column"]] for s in source_rows if s.line==p["source_line"])
    assert original==p["source_value"] and original in p["display_value"]
for full,s in zip(all_full,lookup):
    assert all(full[key]==value for key,value in s.data.items())
assert len(compact["non_hla"])==20 and len(compact["hla"])==4
assert Counter(s.data["endpoint"] for s in lookup)=={"GRAVES_OPHT":12,"E4_GRAVES_OPHT_STRICT":12}
assert all(s.data["match_status"]=="EXACT_GRCH38_REF_ALT" and s.data["qc_status"]=="PASS_GENERAL_QC" for s in lookup)
assert all(p["source_value"] in p["display_value"] for p in provenance)
assert all(sha(ROOT/path)==digest for path,digest in inputs.items() if not path.startswith(rel(DATA)))
dump(DATA/"table_input_manifest.json",{"authoritative_inputs":{p:h for p,h in inputs.items() if not p.startswith(rel(DATA))},"no_statistical_reanalysis":True})
validation={"status":"PASS","no_statistical_reanalysis":True,
    "main_table1_rows":len(table1["rows"]),"main_table1_original_rows_unchanged":True,
    "main_table2_rows":len(table2["rows"]),"main_table2_all_original_row_fields_unchanged":True,
    "hla_lead_rows":len(leads),"hla_same_variant_effect_rows":len(cross),
    "lookup_rows":len(lookup),"lookup_unique_endpoint_variant_pairs":len(test_by_key),
    "non_hla_lookup_rows":20,"hla_lookup_rows":4,"frozen_panel_and_tests_hashes_verified":True,
    "source_field_provenance_records":len(provenance),"all_compact_copied_fields_exact":True,
    "all_original_lookup_columns_exact_in_full_merged_table":True,
    "authoritative_input_hashes_unchanged":True,
    "checks_not_performed":["No statistical model or original-GWAS reanalysis", "No new independent scientific verification", "No Word pagination review by this table builder"]}
dump(REVIEW/"table_build_validation.json",validation)

# Authorized publication-only rounding. The full tables above stay exact.
def rounded(value, specification):
    return format(Decimal(value), specification)


def word_table(ident, title, columns, rows, notes):
    return {"id":ident,"title":title,"columns":[{"key":key,"label":label} for key,label in columns],
            "rows":rows,"notes":notes,"references":[]}


short_trait={"E4_GRAVES_STRICT":"Graves", "E4_HYTHY_AI_STRICT":"AIHT"}
display_note="Publication display only: ORs and approximate normal 95% CI endpoints are rounded to three decimal places, beta and SE to four significant digits, source P and optional adjusted P to three significant digits, and EAF/MAF to three decimal places. Full-precision source values and all original fields are retained in the accompanying machine-readable tables."
word=[]
word.append(word_table("Table S7A", "HLA regional association summary",
    [("trait","Trait"),("raw","Raw rows"),("retained","Retained rows"),("snp","SNPs"),("other","Indel/other"),("p_gw","Source P < 5×10⁻⁸"),("p_suggestive","Source P < 1×10⁻⁵")],
    [{"row_id":s.data["trait_id"],"trait":short_trait[s.data["trait_id"]],"raw":int(s.data["raw_hla_rows"]),"retained":int(s.data["retained_hla_rows"]),"snp":int(s.data["retained_snp_rows"]),"other":int(s.data["retained_indel_or_other_rows"]),"p_gw":int(s.data["source_p_lt_genomewide"]),"p_suggestive":int(s.data["source_p_lt_suggestive"])} for s in summary],
    [hla_note,"The established general QC retains MAF ≥ 0.005. Retained positions span 25,000,462–33,999,677 for both traits; neither boundary has a retained row. Raw source P=0 rows: 0 in each trait. Palindromic retained rows: Graves 8,497; AIHT 8,534. Full regional summary: source_data/table_s7a_region_full.tsv."]))
step_labels={"raw_hla":"Raw HLA records","position_and_nonempty_alleles":"Valid position and nonempty alleles","source_p_finite_gt0_le1":"Finite source P, 0 < P ≤ 1","finite_beta_and_se_gt0":"Finite beta and positive SE","finite_eaf_gt0_lt1":"Finite EAF, 0 < EAF < 1","autosomes":"Autosomal variants","general_maf_min":"MAF ≥ 0.005","deduplicate_variant_id_min_source_p":"Deduplicate variant ID by minimum source P"}
word.append(word_table("Table S7A QC", "Sequential HLA quality control",
    [("trait","Trait"),("step","QC step"),("before","Before"),("after","After"),("removed","Removed")],
    [{"row_id":f"QC{i:02d}","trait":short_trait[s.data["trait_id"]],"step":step_labels[s.data["step"]],"before":int(s.data["rows_before"]),"after":int(s.data["rows_after"]),"removed":int(s.data["removed"])} for i,s in enumerate(flow,1)],
    [hla_note,"Steps retain their original order. All counts are exact. Full QC flow and original step codes: source_data/table_s7a_qc_full.tsv."]))

def word_effect(d, beta="beta", se="se", OR="odds_ratio", low="or_ci95_lower", high="or_ci95_upper", P="source_p", eaf="eaf"):
    return {"beta_se":rounded(d[beta],".4g")+" ("+rounded(d[se],".4g")+")",
            "or_ci":rounded(d[OR],".3f")+" ("+rounded(d[low],".3f")+"–"+rounded(d[high],".3f")+")",
            "source_p":rounded(d[P],".3g"),"frequency":rounded(d[eaf],".3f")+" / "+rounded(d["maf"],".3f")}


word.append(word_table("Table S7B", "Trait-selected HLA regional leads",
    [("selected_for","Selected for"),("variant","Variant (GRCh38)"),("effect_allele","Effect allele"),("beta_se","Beta (SE)"),("or_ci","OR (approximate 95% CI)"),("source_p","Source P"),("frequency","EAF / MAF")],
    [{"row_id":s.data["trait_id"],"selected_for":short_trait[s.data["trait_id"]],"variant":s.data["rsid"]+"\n"+s.data["variant_id"],"effect_allele":s.data["effect_allele"],**word_effect(s.data)} for s in leads],
    [lead_note,"Each trait has one minimum-P lead (one exact tie). The effect allele is ALT; REF is the other allele in the exact variant ID. Full alleles, frequencies, source mlogp and selection metadata: source_data/table_s7b_full.tsv.",ci_note,display_note]))
word.append(word_table("Table S7C", "Same-variant effects at the two HLA regional leads",
    [("selected_for","Selected for"),("evaluated_trait","Evaluated trait"),("variant","Variant (GRCh38)"),("effect_allele","Effect allele"),("beta_se","Beta (SE)"),("or_ci","OR (approximate 95% CI)"),("source_p","Source P"),("frequency","EAF / MAF")],
    [{"row_id":"HLA"+str(i),"selected_for":short_trait[s.data["selected_as_lead_for"]],"evaluated_trait":short_trait[s.data["evaluated_trait"]],"variant":s.data["rsid"]+"\n"+s.data["variant_id"],"effect_allele":s.data["effect_allele"],**word_effect(s.data)} for i,s in enumerate(cross,1)],
    ["All four same-variant marginal estimates are shown. All use exact GRCh38 chromosome, position, REF and ALT matching. No cross-trait effect-difference test was performed.",lead_note,ci_note,display_note,"Full matching, direction, sample counts and source mlogp: source_data/table_s7c_full.tsv."]))
word.append(word_table("Table S8A", "Ophthalmopathy endpoint definitions",
    [("endpoint","Endpoint"),("n_cases","Cases"),("n_controls","Population noncase controls"),("case_rule","Case algorithm"),("control_rule","Verbatim control rule")],
    [{"row_id":s.data["endpoint"],"endpoint":s.data["endpoint"],"n_cases":int(s.data["n_cases"]),"n_controls":int(s.data["n_controls"]),"case_rule":("H7_LACRIMALORBITINOTH condition plus hospital-discharge/cause-of-death ICD-10 registry filters." if s.data["endpoint"]=="GRAVES_OPHT" else "Hospital-discharge/cause-of-death ICD-10 H06.2*E05.9 (raw rule H062&E059)."),"control_rule":s.data["control_definition_verbatim"]} for s in endpoints],
    ["Both endpoints are FinnGen R13, GRCh38. The official labels are Graves ophthalmopathy (GRAVES_OPHT) and Graves ophthalmopathy, strict (E4_GRAVES_OPHT_STRICT)."]+endpoint_notes+
    ["Official definition and control snapshots: "+s.data["official_url"]+"; "+s.data["definition_snapshot"]+", control-rule line "+s.data["definition_line"]+". Counts: "+s.data["count_source"]+", lines "+s.data["count_line"]+"." for s in endpoints]+
    ["Full endpoint metadata and locators: source_data/table_s8a_full.tsv."]))
for endpoint in ("GRAVES_OPHT","E4_GRAVES_OPHT_STRICT"):
    rows=[]
    selected=sorted([s for s in lookup if s.data["endpoint"]==endpoint],key=lambda s:int(s.data["panel_index"]))
    for s in selected:
        d=s.data
        is_hla=int(d["panel_index"])>10
        label=("Graves lead" if d["selection_class"]=="Graves_marginal_HLA_lead" else "AIHT lead") if is_hla else gsu[d["variant_id"]]
        effects=word_effect(d,beta="beta_source_string",se="se_source_string",OR="OR_ALT",low="OR_lower95_normal_approx",high="OR_upper95_normal_approx",P="source_p_string",eaf="eaf_source_string")
        rows.append({"row_id":"P"+d["panel_index"].zfill(2),"group":"HLA" if is_hla else "Non-HLA","variant":d["rsids_source"]+"\n"+label+"\n"+d["variant_id"],"alt":d["alt"],"beta_se":effects["beta_se"],"or_ci":effects["or_ci"],"source_p":effects["source_p"],"optional_p24":rounded(d["p_bonferroni24_optional_descriptive"],".3g")})
    meta=next(s.data for s in endpoints if s.data["endpoint"]==endpoint)
    word.append(word_table("Table S8B1" if endpoint=="GRAVES_OPHT" else "Table S8B2",endpoint+" exploratory associations",
        [("group","Group"),("variant","Variant and selection ID"),("alt","ALT"),("beta_se","Beta (SE)"),("or_ci","OR (approximate 95% CI)"),("source_p","Source P"),("optional_p24","Optional min(1,24×P)")],rows,
        [endpoint+": "+meta["n_cases"]+" cases and "+meta["n_controls"]+" population noncase controls. All 12 frozen variants are shown: ten non-HLA historical discovery representatives and two separately selected HLA regional leads.","All rows have exact GRCh38 REF/ALT matching and pass the established general QC. ALT is the effect allele. Complete QC, matching, missingness, allele frequencies, REF/other alleles, historical discovery context and original precision: source_data/table_s8b_all24_full.tsv."]+lookup_notes[:2]+[ci_note,display_note,"The optional adjustment uses the complete family of 24 lookups across both endpoints. It is a descriptive threshold after discovery selection, not independent replication.",freeze_note]))
assert [len(p["rows"]) for p in word]==[2,16,2,4,2,12,12]
dump(TABLES/"word_supplement_tables.json",word)
word_index=[]
for i,p in enumerate(word):
    name="word_"+("s7a_region","s7a_qc","s7b","s7c","s8a","s8b_graves_opht","s8b_strict")[i]
    tsv(TABLES/f"{name}.tsv",["row_id"]+[c["key"] for c in p["columns"]],p["rows"])
    word_index.append({"index":i,"id":p["id"],"title":p["title"],"row_count":len(p["rows"]),"column_keys":[c["key"] for c in p["columns"]],"tsv":rel(TABLES/f"{name}.tsv")})
dump(TABLES/"word_supplement_table_index.json",word_index)
validation["publication_display_views"]=7
validation["publication_display_lookup_rows"]=24
validation["publication_display_rounding_authorized"]=True
dump(REVIEW/"table_build_validation.json",validation)
print(json.dumps(validation,indent=2))
