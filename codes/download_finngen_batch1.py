#!/usr/bin/env python3
"""
FinnGen Data Downloader - Batch 1 (Core GWAS & Multiome)
Autoimmune Thyroid Disease Project
"""

import os
import sys
import csv
import json
import hashlib
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DATA_DIR = os.path.join(BASE_DIR, "raw_data")
MANIFEST_DIR = os.path.join(RAW_DATA_DIR, "manifests")
LOCAL_MANIFEST_PATH = os.path.join(MANIFEST_DIR, "local_data_manifest.tsv")

MANIFEST_HEADERS = [
    "source", "release", "dataset_type", "trait_or_celltype", "filename",
    "url_or_gcs_path", "file_size", "checksum", "genome_build", "effect_allele",
    "gene_annotation", "download_date", "citation", "license"
]

def get_md5(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def download_file(url, local_path, expected_size=None, expected_md5=None, max_retries=5):
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    
    # Check if already downloaded and valid
    if os.path.exists(local_path):
        cur_size = os.path.getsize(local_path)
        if expected_size is not None:
            if cur_size == expected_size:
                if expected_md5:
                    cur_md5 = get_md5(local_path)
                    if cur_md5 == expected_md5:
                        print(f"[SKIP] Already exists and verified: {os.path.basename(local_path)}")
                        return True, local_path, cur_size, cur_md5
                else:
                    print(f"[SKIP] Already exists with correct size: {os.path.basename(local_path)}")
                    return True, local_path, cur_size, get_md5(local_path)
        else:
            if cur_size > 1024:
                print(f"[SKIP] Already exists: {os.path.basename(local_path)} ({cur_size / (1024*1024):.2f} MB)")
                return True, local_path, cur_size, get_md5(local_path)

    for attempt in range(1, max_retries + 1):
        try:
            print(f"[DOWNLOADING] {os.path.basename(local_path)} (Attempt {attempt}/{max_retries})")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=60) as resp, open(local_path, "wb") as out_file:
                chunk_size = 1024 * 1024 * 2 # 2MB
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    out_file.write(chunk)
            
            actual_size = os.path.getsize(local_path)
            actual_md5 = get_md5(local_path)
            
            if expected_size is not None and actual_size != expected_size:
                print(f"[WARN] Size mismatch for {local_path}: expected {expected_size}, got {actual_size}")
                time.sleep(2)
                continue
            if expected_md5 and actual_md5 != expected_md5:
                print(f"[WARN] MD5 mismatch for {local_path}: expected {expected_md5}, got {actual_md5}")
                time.sleep(2)
                continue
                
            print(f"[SUCCESS] Downloaded {os.path.basename(local_path)} ({actual_size / (1024*1024):.2f} MB)")
            return True, local_path, actual_size, actual_md5
        except Exception as e:
            print(f"[ERROR] Attempt {attempt} failed for {url}: {e}")
            time.sleep(3)
            
    return False, local_path, 0, ""

def append_to_local_manifest(entries):
    file_exists = os.path.exists(LOCAL_MANIFEST_PATH)
    existing_files = set()
    if file_exists:
        with open(LOCAL_MANIFEST_PATH, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                existing_files.add(row["filename"])
                
    with open(LOCAL_MANIFEST_PATH, "a" if file_exists else "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_HEADERS, delimiter="\t")
        if not file_exists:
            writer.writeheader()
        for entry in entries:
            if entry["filename"] not in existing_files:
                writer.writerow(entry)
                existing_files.add(entry["filename"])

def build_batch1_plan():
    tasks = []
    
    # -------------------------------------------------------------
    # 1. FinnGen R13 GWAS Summary Stats
    # -------------------------------------------------------------
    r13_endpoints = [
        ("E4_GRAVES_STRICT", "Graves_strict"),
        ("E4_HYTHY_AI_STRICT", "AIHT_strict"),
        ("AUTOIMMUNE_HYPERTHYROIDISM", "Autoimmune_hyperthyroidism"),
        ("E4_THYTOXGOITDIF", "Thyrotoxicosis_diffuse_goitre"),
        ("GRAVES_OPHT", "Graves_ophthalmopathy"),
        ("E4_GRAVES_OPHT_STRICT", "Graves_ophthalmopathy_strict"),
        ("AUTOIMMUNE", "Autoimmune_all"),
        ("AUTOIMMUNE_NONTHYROID", "Autoimmune_nonthyroid")
    ]
    
    for ep, label in r13_endpoints:
        # GWAS .gz
        url_gz = f"https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{ep}.gz"
        local_gz = os.path.join(RAW_DATA_DIR, "finngen_r13", "summary_stats", f"finngen_R13_{ep}.gz")
        tasks.append({
            "source": "FinnGen", "release": "R13", "dataset_type": "GWAS_SumStats",
            "trait_or_celltype": ep, "filename": f"finngen_R13_{ep}.gz",
            "url": url_gz, "local_path": local_gz, "expected_size": None, "expected_md5": None,
            "genome_build": "GRCh38", "effect_allele": "ALT", "gene_annotation": "GENCODE_v38/Ensembl",
            "citation": "Kurki et al. Nature 2023; FinnGen R13 2026", "license": "FinnGen Public Access"
        })
        # GWAS .tbi
        url_tbi = f"https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{ep}.gz.tbi"
        local_tbi = os.path.join(RAW_DATA_DIR, "finngen_r13", "summary_stats", f"finngen_R13_{ep}.gz.tbi")
        tasks.append({
            "source": "FinnGen", "release": "R13", "dataset_type": "GWAS_Index",
            "trait_or_celltype": ep, "filename": f"finngen_R13_{ep}.gz.tbi",
            "url": url_tbi, "local_path": local_tbi, "expected_size": None, "expected_md5": None,
            "genome_build": "GRCh38", "effect_allele": "ALT", "gene_annotation": "GENCODE_v38/Ensembl",
            "citation": "Kurki et al. Nature 2023; FinnGen R13 2026", "license": "FinnGen Public Access"
        })
        
        # SuSiE fine-mapping
        for ext in ["SUSIE.cred.bgz", "SUSIE.cred_99.bgz", "SUSIE.snp.bgz", "SUSIE.snp.bgz.tbi"]:
            url_fm = f"https://storage.googleapis.com/finngen-public-data-r13/finemap/full/susie/finngen_R13_{ep}.{ext}"
            local_fm = os.path.join(RAW_DATA_DIR, "finngen_r13", "finemap", "susie", f"finngen_R13_{ep}.{ext}")
            tasks.append({
                "source": "FinnGen", "release": "R13", "dataset_type": "GWAS_SuSiE_FineMapping",
                "trait_or_celltype": ep, "filename": f"finngen_R13_{ep}.{ext}",
                "url": url_fm, "local_path": local_fm, "expected_size": None, "expected_md5": None,
                "genome_build": "GRCh38", "effect_allele": "ALT", "gene_annotation": "GENCODE_v38/Ensembl",
                "citation": "Kurki et al. Nature 2023; FinnGen R13 2026", "license": "FinnGen Public Access"
            })
            
        # FINEMAP fine-mapping
        for ext in ["FINEMAP.config.bgz", "FINEMAP.region.bgz", "FINEMAP.snp.bgz", "FINEMAP.snp.bgz.tbi"]:
            url_fm = f"https://storage.googleapis.com/finngen-public-data-r13/finemap/full/finemap/finngen_R13_{ep}.{ext}"
            local_fm = os.path.join(RAW_DATA_DIR, "finngen_r13", "finemap", "finemap", f"finngen_R13_{ep}.{ext}")
            tasks.append({
                "source": "FinnGen", "release": "R13", "dataset_type": "GWAS_FINEMAP_FineMapping",
                "trait_or_celltype": ep, "filename": f"finngen_R13_{ep}.{ext}",
                "url": url_fm, "local_path": local_fm, "expected_size": None, "expected_md5": None,
                "genome_build": "GRCh38", "effect_allele": "ALT", "gene_annotation": "GENCODE_v38/Ensembl",
                "citation": "Kurki et al. Nature 2023; FinnGen R13 2026", "license": "FinnGen Public Access"
            })

    # -------------------------------------------------------------
    # 2. FinnGen Multiome Files from Manifest
    # -------------------------------------------------------------
    multiome_manifest = os.path.join(MANIFEST_DIR, "finngen_multiome_v1.manifest.tsv")
    with open(multiome_manifest, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rel_path = row["path"]
            product = row["product"]
            analysis = row["analysis"]
            sz = int(row["size_bytes"])
            md5 = row["md5"]
            celltype = row.get("celltype", "")
            
            # Decide if included in Batch 1 core download
            # Included in Batch 1:
            # - atac_peaks (consensus_peaks only, skip huge BigWigs in Batch 1)
            # - cascade (all 36 tables)
            # - coloc (all 4 tables)
            # - peak_gene_links (all 33 analysis set files)
            # - eQTL/cis (all 33 cis_acat files)
            # - caQTL/cis (all 33 cis_acat files)
            # - eQTL/finemapping (*.SUSIE.cred*.tsv.gz - credible sets summary)
            # - caQTL/finemapping (*.SUSIE.cred*.tsv.gz - credible sets summary)
            # - cis_nominal (variant-level QTL sumstats for colocalization/MR)
            
            include = False
            dtype = ""
            
            if product == "atac_peaks" and analysis == "consensus_peaks":
                include = True
                dtype = "Multiome_Consensus_Peaks"
            elif product == "cascade":
                include = True
                dtype = "Multiome_CASCADE_Categorization"
            elif product == "coloc":
                include = True
                dtype = "Multiome_Coloc_Baseline"
            elif product == "peak_gene_links" and analysis == "":
                include = True
                dtype = "Multiome_Peak_Gene_Links"
            elif product == "cis" and analysis == "cis_acat":
                include = True
                dtype = f"Multiome_{row.get('modality', '')}_ACAT"
            elif product == "finemapping" and ("cred" in analysis):
                include = True
                dtype = f"Multiome_{row.get('modality', '')}_SuSiE_CredSet"
            elif product == "cis_nominal":
                include = True
                dtype = f"Multiome_{row.get('modality', '')}_Nominal_SumStats"
                
            if include:
                url = f"https://storage.googleapis.com/finngen-public-data-multiome/{rel_path}"
                local_path = os.path.join(RAW_DATA_DIR, "finngen_multiome", rel_path)
                tasks.append({
                    "source": "FinnGen", "release": "Multiome_v1_2026", "dataset_type": dtype,
                    "trait_or_celltype": celltype, "filename": os.path.basename(rel_path),
                    "url": url, "local_path": local_path, "expected_size": sz, "expected_md5": md5,
                    "genome_build": "GRCh38", "effect_allele": "allele2 (ALT)", "gene_annotation": "GENCODE_v38/Ensembl",
                    "citation": "Kanai et al. medRxiv 2025 (doi:10.1101/2025.11.25.25340489)", "license": "FinnGen Public Access"
                })
                
    return tasks

if __name__ == "__main__":
    tasks = build_batch1_plan()
    print(f"Total Batch 1 tasks planned: {len(tasks)}")
    
    # Separate small/fast metadata and large files
    fast_tasks = [t for t in tasks if t["dataset_type"] != "Multiome_eQTL_Nominal_SumStats" and t["dataset_type"] != "Multiome_caQTL_Nominal_SumStats" and t["dataset_type"] != "GWAS_SumStats"]
    large_tasks = [t for t in tasks if t not in fast_tasks]
    
    print(f"  Fast/Medium tasks (GWAS fine-mapping, Multiome CASCADE, links, coloc, peaks, ACAT, credible sets): {len(fast_tasks)}")
    print(f"  Large tasks (GWAS full sumstats, Multiome cis_nominal): {len(large_tasks)}")
    
    today_str = time.strftime("%Y-%m-%d")
    
    # 1. Execute fast tasks first (multi-threaded)
    print("\n>>> Phase 1: Downloading Core Summaries, Fine-mapping, CASCADE, Links, Peaks, Coloc...")
    manifest_entries = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        future_to_task = {executor.submit(download_file, t["url"], t["local_path"], t["expected_size"], t["expected_md5"]): t for t in fast_tasks}
        for future in as_completed(future_to_task):
            t = future_to_task[future]
            success, path, sz, md5_val = future.result()
            if success:
                manifest_entries.append({
                    "source": t["source"], "release": t["release"], "dataset_type": t["dataset_type"],
                    "trait_or_celltype": t["trait_or_celltype"], "filename": t["filename"],
                    "url_or_gcs_path": t["url"], "file_size": sz, "checksum": md5_val,
                    "genome_build": t["genome_build"], "effect_allele": t["effect_allele"],
                    "gene_annotation": t["gene_annotation"], "download_date": today_str,
                    "citation": t["citation"], "license": t["license"]
                })
                
    append_to_local_manifest(manifest_entries)
    print(f"Phase 1 complete. Manifest updated with {len(manifest_entries)} entries.")
    
    # 2. Execute large tasks (GWAS sumstats first, then multiome nominals)
    print("\n>>> Phase 2: Downloading GWAS Full Summary Statistics & Multiome Nominal SumStats...")
    gwas_tasks = [t for t in large_tasks if t["dataset_type"] == "GWAS_SumStats"]
    multiome_nom_tasks = [t for t in large_tasks if "Nominal" in t["dataset_type"]]
    
    # Download GWAS sumstats (8 files, ~6GB total)
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_task = {executor.submit(download_file, t["url"], t["local_path"], t["expected_size"], t["expected_md5"]): t for t in gwas_tasks}
        for future in as_completed(future_to_task):
            t = future_to_task[future]
            success, path, sz, md5_val = future.result()
            if success:
                append_to_local_manifest([{
                    "source": t["source"], "release": t["release"], "dataset_type": t["dataset_type"],
                    "trait_or_celltype": t["trait_or_celltype"], "filename": t["filename"],
                    "url_or_gcs_path": t["url"], "file_size": sz, "checksum": md5_val,
                    "genome_build": t["genome_build"], "effect_allele": t["effect_allele"],
                    "gene_annotation": t["gene_annotation"], "download_date": today_str,
                    "citation": t["citation"], "license": t["license"]
                }])
    print("Phase 2 (GWAS SumStats) complete.")
    
    # Download Multiome nominal sumstats
    print(f"\n>>> Phase 3: Downloading Multiome Nominal SumStats ({len(multiome_nom_tasks)} files)...")
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_task = {executor.submit(download_file, t["url"], t["local_path"], t["expected_size"], t["expected_md5"]): t for t in multiome_nom_tasks}
        for future in as_completed(future_to_task):
            t = future_to_task[future]
            success, path, sz, md5_val = future.result()
            if success:
                append_to_local_manifest([{
                    "source": t["source"], "release": t["release"], "dataset_type": t["dataset_type"],
                    "trait_or_celltype": t["trait_or_celltype"], "filename": t["filename"],
                    "url_or_gcs_path": t["url"], "file_size": sz, "checksum": md5_val,
                    "genome_build": t["genome_build"], "effect_allele": t["effect_allele"],
                    "gene_annotation": t["gene_annotation"], "download_date": today_str,
                    "citation": t["citation"], "license": t["license"]
                }])
    print("Batch 1 Download and local manifest generation completely finished!")
