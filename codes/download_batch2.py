#!/usr/bin/env python3
"""
Batch 2 External Datasets Downloader
Autoimmune Thyroid Disease Project
"""

import os
import sys
import csv
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

def is_file_valid(file_path):
    if not os.path.exists(file_path):
        return False
    if os.path.getsize(file_path) == 0:
        return False
    if file_path.endswith(".tar.gz") or file_path.endswith(".tgz"):
        try:
            import tarfile
            with tarfile.open(file_path, "r:gz") as f:
                _ = f.getmembers()
            return True
        except Exception:
            return False
    elif file_path.endswith(".tar"):
        try:
            import tarfile
            with tarfile.open(file_path, "r:") as f:
                _ = f.getmembers()
            return True
        except Exception:
            return False
    elif file_path.endswith(".gz"):
        try:
            import gzip
            with gzip.open(file_path, "rb") as f:
                while f.read(1024 * 1024 * 5):
                    pass
            return True
        except Exception:
            return False
    return True

def download_file(task, max_retries=10):
    url = task["url"]
    local_path = task["local_path"]
    part_path = local_path + ".part"
    
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    
    if is_file_valid(local_path):
        cur_size = os.path.getsize(local_path)
        cur_md5 = get_md5(local_path)
        print(f"[SKIP] Already exists and verified valid: {os.path.basename(local_path)} ({cur_size / (1024*1024):.2f} MB)")
        return True, task, cur_size, cur_md5

    for attempt in range(1, max_retries + 1):
        try:
            print(f"[DOWNLOADING] {os.path.basename(local_path)} (Attempt {attempt}/{max_retries})")
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=180) as resp, open(part_path, "wb") as out_file:
                content_length = resp.headers.get("Content-Length")
                expected_len = int(content_length) if content_length and content_length.isdigit() else None
                
                chunk_size = 1024 * 1024 * 2
                downloaded = 0
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    out_file.write(chunk)
                    downloaded += len(chunk)
                    
            if expected_len is not None and downloaded != expected_len:
                print(f"[WARN] Length mismatch for {part_path}: expected {expected_len}, got {downloaded}")
                time.sleep(2)
                continue
                
            if os.path.exists(part_path) and is_file_valid(part_path):
                if os.path.exists(local_path):
                    os.remove(local_path)
                os.rename(part_path, local_path)
                actual_size = os.path.getsize(local_path)
                actual_md5 = get_md5(local_path)
                print(f"[SUCCESS] Downloaded & Verified: {os.path.basename(local_path)} ({actual_size / (1024*1024):.2f} MB)")
                return True, task, actual_size, actual_md5
            else:
                print(f"[WARN] Integrity check failed for {part_path}, retrying...")
                time.sleep(2)
        except Exception as e:
            print(f"[ERROR] Attempt {attempt} failed for {url}: {e}")
            time.sleep(3)
            
    return False, task, 0, ""

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

def build_batch2_plan():
    tasks = []
    
    # -------------------------------------------------------------
    # 1. GEO GSE182416 (正常甲状腺单细胞参考)
    # -------------------------------------------------------------
    tasks.append({
        "source": "NCBI_GEO", "release": "GSE182416", "dataset_type": "scRNA_CountMatrix",
        "trait_or_celltype": "Normal_Thyroid", "filename": "GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182416/suppl/GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz",
        "local_path": os.path.join(RAW_DATA_DIR, "single_cell", "GSE182416", "GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz"),
        "expected_size": None, "expected_md5": None,
        "genome_build": "GRCh38", "effect_allele": "NA", "gene_annotation": "GENCODE/Ensembl",
        "citation": "GSE182416 (Thyroid Single-cell Atlas)", "license": "NCBI Open Data"
    })
    tasks.append({
        "source": "NCBI_GEO", "release": "GSE182416", "dataset_type": "scRNA_Metadata",
        "trait_or_celltype": "Normal_Thyroid", "filename": "GSE182416_Thyroid_normal_7samples_metadata.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182416/suppl/GSE182416_Thyroid_normal_7samples_metadata.txt.gz",
        "local_path": os.path.join(RAW_DATA_DIR, "single_cell", "GSE182416", "GSE182416_Thyroid_normal_7samples_metadata.txt.gz"),
        "expected_size": None, "expected_md5": None,
        "genome_build": "GRCh38", "effect_allele": "NA", "gene_annotation": "GENCODE/Ensembl",
        "citation": "GSE182416 (Thyroid Single-cell Atlas)", "license": "NCBI Open Data"
    })

    # -------------------------------------------------------------
    # 2. GEO GSE248205 (AITD 空间转录组 Visium)
    # -------------------------------------------------------------
    tasks.append({
        "source": "NCBI_GEO", "release": "GSE248205", "dataset_type": "Spatial_Visium_Processed",
        "trait_or_celltype": "AITD_Spatial", "filename": "GSE248205_Processed_data.tar.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248205/suppl/GSE248205_Processed_data.tar.gz",
        "local_path": os.path.join(RAW_DATA_DIR, "spatial", "GSE248205", "GSE248205_Processed_data.tar.gz"),
        "expected_size": None, "expected_md5": None,
        "genome_build": "GRCh38", "effect_allele": "NA", "gene_annotation": "10x SpaceRanger",
        "citation": "GSE248205 (AITD Spatial Transcriptomics)", "license": "NCBI Open Data"
    })

    # -------------------------------------------------------------
    # 3. GEO GSE285196 (Graves B 细胞 scRNA/BCR)
    # -------------------------------------------------------------
    tasks.append({
        "source": "NCBI_GEO", "release": "GSE285196", "dataset_type": "scRNA_BCR_Matrix",
        "trait_or_celltype": "Graves_B_cells", "filename": "GSE285196_RAW.tar",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285196/suppl/GSE285196_RAW.tar",
        "local_path": os.path.join(RAW_DATA_DIR, "single_cell", "GSE285196", "GSE285196_RAW.tar"),
        "expected_size": None, "expected_md5": None,
        "genome_build": "GRCh38", "effect_allele": "NA", "gene_annotation": "10x Genomics CellRanger",
        "citation": "GSE285196 (Graves B Cell Transcriptomics)", "license": "NCBI Open Data"
    })

    # -------------------------------------------------------------
    # 4. GEO GSE218743 (Hashimoto CD45+ 免疫单细胞)
    # -------------------------------------------------------------
    tasks.append({
        "source": "NCBI_GEO", "release": "GSE218743", "dataset_type": "scRNA_CountMatrix",
        "trait_or_celltype": "Hashimoto_CD45_plus", "filename": "GSE218743_irae_ht_integrated_counts.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218743/suppl/GSE218743_irae_ht_integrated_counts.csv.gz",
        "local_path": os.path.join(RAW_DATA_DIR, "single_cell", "GSE218743", "GSE218743_irae_ht_integrated_counts.csv.gz"),
        "expected_size": None, "expected_md5": None,
        "genome_build": "GRCh38", "effect_allele": "NA", "gene_annotation": "GENCODE/Ensembl",
        "citation": "GSE218743 (Hashimoto CD45+ scRNA-seq)", "license": "NCBI Open Data"
    })
    tasks.append({
        "source": "NCBI_GEO", "release": "GSE218743", "dataset_type": "scRNA_Metadata",
        "trait_or_celltype": "Hashimoto_CD45_plus", "filename": "GSE218743_irae_ht_integrated_metadata.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218743/suppl/GSE218743_irae_ht_integrated_metadata.csv.gz",
        "local_path": os.path.join(RAW_DATA_DIR, "single_cell", "GSE218743", "GSE218743_irae_ht_integrated_metadata.csv.gz"),
        "expected_size": None, "expected_md5": None,
        "genome_build": "GRCh38", "effect_allele": "NA", "gene_annotation": "GENCODE/Ensembl",
        "citation": "GSE218743 (Hashimoto CD45+ scRNA-seq)", "license": "NCBI Open Data"
    })

    return tasks

if __name__ == "__main__":
    tasks = build_batch2_plan()
    print(f"Total Batch 2 tasks planned: {len(tasks)}")
    
    today_str = time.strftime("%Y-%m-%d")
    manifest_entries = []
    
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_task = {executor.submit(download_file, t): t for t in tasks}
        for future in as_completed(future_to_task):
            t = future_to_task[future]
            success, task_info, sz, md5_val = future.result()
            if success:
                manifest_entries.append({
                    "source": task_info["source"], "release": task_info["release"], "dataset_type": task_info["dataset_type"],
                    "trait_or_celltype": task_info["trait_or_celltype"], "filename": task_info["filename"],
                    "url_or_gcs_path": task_info["url"], "file_size": sz, "checksum": md5_val,
                    "genome_build": task_info["genome_build"], "effect_allele": task_info["effect_allele"],
                    "gene_annotation": task_info["gene_annotation"], "download_date": today_str,
                    "citation": task_info["citation"], "license": task_info["license"]
                })
                
    append_to_local_manifest(manifest_entries)
    print(f"Batch 2 download complete! Manifest updated with {len(manifest_entries)} entries.")
