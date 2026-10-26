#!/usr/bin/env python3
"""Compile the source-curated T15 registry and conservative cohort overlap gate.

Standard library only. Does not download GWAS statistics or run replication.
Missing participant overlap counts stay missing, never zero.
"""
import csv
import hashlib
import itertools
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outcome_data/t15"
CONFIG = ROOT / "config/t15_external_cohorts.json"
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: ";".join(row[k]) if isinstance(row.get(k), list)
                        else row.get(k) for k in fields})


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def xlsx_cells(path, sheet_name):
    """Read exact source cells without changing the publisher workbook."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            strings = ["".join(t.text or "" for t in si.iterfind(".//s:t", NS))
                       for si in root.findall("s:si", NS)]
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheet = next(s for s in wb.findall("s:sheets/s:sheet", NS) if s.attrib["name"] == sheet_name)
        rid = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(r.attrib["Target"] for r in rels if r.attrib["Id"] == rid)
        member = target.lstrip("/") if target.startswith("/") else str(PurePosixPath("xl") / target)
        data = {}
        for cell in ET.fromstring(z.read(member)).findall(".//s:sheetData/s:row/s:c", NS):
            value = cell.find("s:v", NS)
            text = value.text if value is not None else None
            if cell.attrib.get("t") == "s" and text is not None:
                text = strings[int(text)]
            elif cell.attrib.get("t") == "inlineStr":
                text = "".join(t.text or "" for t in cell.findall(".//s:t", NS))
            data[cell.attrib["r"]] = text
        return data


def relationship(a, b):
    if a["dataset_id"] == b["dataset_id"]:
        return "SELF", [], "same registered dataset"
    if not (a["cohort_membership_complete"] and b["cohort_membership_complete"]):
        return "UNCONFIRMED", [], "exact source/subset membership unresolved"
    common = sorted(set(a["cohorts"]) & set(b["cohorts"]))
    if common:
        return "PARTIAL_OVERLAP", common, "shared cohort sampling frame; exact participant overlap unknown"
    if "23ANDME" in a["cohorts"] + b["cohorts"]:
        return "UNCONFIRMED", [], "23andMe cross-enrollment not ruled out"
    return "NO_SHARED_COHORT_IDENTIFIED", [], "cohort-level inference only; no participant-level linkage performed"


def main():
    cfg = json.loads(CONFIG.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    sources = {s["source_id"]: s for s in cfg["sources"]}
    cells = xlsx_cells(ROOT / sources["S02"]["local_path"], "ST1")
    canonical = ["ALL_OF_US", "BBJ", "CROHT", "DISCOVEHR", "ESTBB",
                 "GODARTS_SHARE", "BMC_LATVIA", "HUNT", "MGI", "UKB"]
    cohorts = []
    for col, cohort_id in zip("CDEFGHIJKL", canonical):
        n_cases, n_controls = int(cells[col + "6"]), int(cells[col + "7"])
        if cohort_id in ("ALL_OF_US", "UKB"):
            pattern = r"EUR(?: \(PCA-based\))?:\s*(\d+)"
            n_cases_eur = int(re.search(pattern, cells[col + "11"]).group(1))
            n_controls_eur = int(re.search(pattern, cells[col + "12"]).group(1))
        elif cohort_id == "BBJ":
            n_cases_eur, n_controls_eur = 0, 0
        else:
            n_cases_eur, n_controls_eur = n_cases, n_controls
        cohorts.append({
            "cohort_id": cohort_id, "source_name": cells[col + "3"],
            "ancestry_as_reported": cells[col + "4"], "country": cells[col + "5"],
            "n_cases": n_cases, "n_controls": n_controls,
            "n_cases_eur": n_cases_eur, "n_controls_eur": n_controls_eur,
            "component_build_as_reported": cells[col + "17"],
            "meta_output_build": "GRCh37",
            "source_id": "S02", "source_cells": f"ST1!{col}3:{col}7;{col}11:{col}12;{col}17",
        })
    all_cohorts = [r["cohort_id"] for r in cohorts]
    eur_cohorts = [c for c in all_cohorts if c != "BBJ"]
    aliases = {"AUTO_THYRO_ALL": all_cohorts, "AUTO_THYRO_EUR": eur_cohorts,
               "AUTO_THYRO_EUR_NON_UKB": [c for c in eur_cohorts if c != "UKB"]}
    rows = cfg["datasets"]
    for r in rows:
        r["cohorts"] = sorted(c for c0 in r["cohorts"] for c in aliases.get(c0, [c0]))
        complete = r["cohort_membership_complete"]
        r["overlap_vs_finngen"] = ("UNCONFIRMED" if not complete else
                                   "PARTIAL_OVERLAP" if "FINNGEN" in r["cohorts"] else "NO_SHARED_COHORT_IDENTIFIED")
        r["overlap_vs_ukb"] = ("UNCONFIRMED" if not complete else
                               "PARTIAL_OVERLAP" if "UKB" in r["cohorts"] else "NO_SHARED_COHORT_IDENTIFIED")
        r["strict_independence_eligible"] = complete and not (
            set(r["cohorts"]) & set(cfg["strict_independence_excluded_cohorts"]))
        r["metadata_complete_for_acquisition"] = all(
            r[k] is not None for k in ("n_cases", "n_controls")) and r["genome_build"] in ("GRCh37", "GRCh38")
        r["replication_executed"] = False
        r["source_urls"] = [sources[s]["url"] for s in r["source_ids"]]
    by_id = {r["dataset_id"]: r for r in rows}
    edges = []
    matrix = []
    labels = {"SELF": "自身", "PARTIAL_OVERLAP": "部分重叠",
              "UNCONFIRMED": "不可确认", "NO_SHARED_COHORT_IDENTIFIED": "完全独立_仅队列组成层面"}
    for a in rows:
        matrix.append({"dataset_id": a["dataset_id"], **{
            b["dataset_id"]: relationship(a, b)[0] for b in rows}})
    for a, b in itertools.combinations(rows, 2):
        status, shared, basis = relationship(a, b)
        edges.append({"dataset_a": a["dataset_id"], "dataset_b": b["dataset_id"],
                      "overlap_class": status, "overlap_class_zh": labels[status],
                      "shared_cohorts": shared, "shared_cases_n": None, "shared_controls_n": None,
                      "basis": basis, "source_ids": sorted(set(a["source_ids"] + b["source_ids"]))})
    traits = {r["trait_id"]: r for r in read_tsv(ROOT / "config/traits.tsv")}
    gd_catalog = json.loads((OUT / "sources/GCST90018847.json").read_text())
    ht_catalog = json.loads((OUT / "sources/GCST90018855.json").read_text())
    bbj_catalog = json.loads((OUT / "sources/GCST90018627.json").read_text())
    def eur_countries(d):
        return {c["countryName"] for a in d["ancestries"]
                if any(g["ancestralGroup"] == "European" for g in a["ancestralGroups"])
                for c in a["countryOfRecruitment"]}
    checks = {
        "dataset_ids_unique": len(by_id) == len(rows),
        "source_ids_resolve": all(s in sources for r in rows for s in r["source_ids"]),
        "all_metadata_fields_explicit": all(k in r for r in rows for k in
            ("ancestry", "n_cases", "n_controls", "phenotype", "genome_build", "data_availability")),
        "thyroidomics_ten_source_cohorts": len(cohorts) == 10 and len(set(all_cohorts)) == 10,
        "thyroidomics_case_control_sums": sum(r["n_cases"] for r in cohorts) == 48694 and sum(r["n_controls"] for r in cohorts) == 1044134,
        "thyroidomics_eur_counts": sum(c["n_cases_eur"] for c in cohorts) == by_id["THYRO2026_EUR"]["n_cases"] and sum(c["n_controls_eur"] for c in cohorts) == by_id["THYRO2026_EUR"]["n_controls"],
        "thyroidomics_ukb_present_finngen_absent": "UKB" in all_cohorts and "FINNGEN" not in all_cohorts,
        "sakaue_eur_finland_uk_confirmed": all({"Finland", "U.K."} <= eur_countries(d) for d in [gd_catalog, ht_catalog]),
        "bbj2021_source_is_bbj_only": bbj_catalog["cohort"] == "BBJ",
        "frozen_discovery_counts": all(int(traits[t][k]) == by_id[i][k] for t,i in
            [("E4_HYTHY_AI_STRICT", "FG13_AIHT"), ("E4_GRAVES_STRICT", "FG13_GD")] for k in ["n_cases", "n_controls"]),
        "matrix_complete_symmetric": len(edges) == len(rows)*(len(rows)-1)//2 and all(
            matrix[i][b["dataset_id"]] == matrix[j][a["dataset_id"]]
            for i,a in enumerate(rows) for j,b in enumerate(rows)),
        "no_shared_participant_counts_invented": all(e["shared_cases_n"] is None and e["shared_controls_n"] is None for e in edges),
        "no_finngen_or_ukb_in_independent_designation": all(not r["strict_independence_eligible"] for r in rows if set(r["cohorts"]) & {"FINNGEN", "UKB"}),
        "selected_graves_plan_independent": by_id["BBJ_2019_GD"]["strict_independence_eligible"],
        "bbj_releases_not_two_independent_replications": relationship(by_id["BBJ_2019_GD"], by_id["BBJ_2021_GD"])[0] == "PARTIAL_OVERLAP",
        "requested_split_not_claimed_available": not by_id["THYRO2026_NON_UKB_EUR_REQUEST"]["summary_statistics_downloaded"] and not by_id["THYRO2026_NON_UKB_EUR_REQUEST"]["cohort_membership_complete"],
        "t12_t13_t14_holds_preserved": all(read_tsv(ROOT / p)[0]["status"] == s for p,s in
            [("outcome_data/t12/t12_final_decision.tsv", "HOLD"), ("outcome_data/t13/t13_decision.tsv", "PASS_WITH_HOLD"), ("outcome_data/t14/t14_decision.tsv", "PASS_WITH_HOLD")]),
        "t16_t17_not_executed": not any(r["replication_executed"] for r in rows),
    }
    tsv(OUT / "t15_external_cohort_registry.tsv", rows)
    tsv(OUT / "t15_thyroidomics_component_cohorts.tsv", cohorts)
    tsv(OUT / "t15_overlap_matrix.tsv", matrix)
    tsv(OUT / "t15_overlap_pairs.tsv", edges)
    tsv(OUT / "t15_sources.tsv", list(sources.values()))
    tsv(OUT / "t15_validation.tsv", [{"check":k,"passed":v} for k,v in checks.items()])
    provenance = []
    paths = [CONFIG, Path(__file__), ROOT / "config/traits.tsv", ROOT / "config/analysis.yml"]
    paths += [ROOT / s["local_path"] for s in sources.values() if s["local_path"]]
    paths += [ROOT / p for p in ["outcome_data/t12/t12_final_decision.tsv", "outcome_data/t13/t13_decision.tsv", "outcome_data/t14/t14_decision.tsv"]]
    for p in sorted(set(paths)):
        provenance.append({"path":str(p.relative_to(ROOT)),"bytes":p.stat().st_size,"sha256":sha256(p)})
    tsv(OUT / "t15_provenance.tsv", provenance)
    result = {
        "gate":"T15_external_cohort_and_overlap_audit",
        "status":"PASS_WITH_HOLD" if all(checks.values()) else "FAIL",
        "decision":"PASS_T15_REGISTRY_HOLD_INDEPENDENT_AIHT_AND_METADATA" if all(checks.values()) else "FAIL_T15_VALIDATION",
        "datasets":len(rows),"pairwise_relationships":len(edges),"checks_passed":sum(checks.values()),"checks_total":len(checks),
        "graves_plan":"BBJ_2019_GD; cross-ancestry; T17 not executed",
        "aiht_plan":"Request non-FinnGen/non-UKB component or split; ThyroidOmics public metas comparison-only under frozen gate",
        "strict_independence_excludes":"FINNGEN;UKB",
        "exact_person_overlap":"UNKNOWN; no individual-level identity linkage",
        "completed_at_utc":datetime.now(timezone.utc).isoformat(),
    }
    tsv(OUT / "t15_decision.tsv", [result])
    print(json.dumps({**result,"pair_classes":dict(Counter(e["overlap_class"] for e in edges))}, ensure_ascii=False, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
