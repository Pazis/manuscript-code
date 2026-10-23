"""Fetch, select and download cyanobacterial genomes from NCBI.

Run in the directory where the files should go:
  python genome_fetcher/main.py {fetch,filter,download,organize}
"""
import argparse
import shutil
import sys
import time
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

TAXON_ID = 1117  # Cyanobacteriota
OUTPUT_RAW = "cyanobacteria_raw.csv"
OUTPUT_FILTERED = "cyanobacteria_selected.xlsx"
ACCESSION_LIST = "accession_list.txt"
ZIP_FILE = "genomes_dataset.zip"
FINAL_DATASET_DIR = "Cyano_Benchmark_Dataset"

API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def get_lineage_batch(taxids):
    """Family and order per TaxID from NCBI Taxonomy. Failed batches are reported and left out."""
    results = {}
    batch_size = 100
    print(f"Fetching taxonomy for {len(taxids)} unique TaxIDs...")

    for i in tqdm(range(0, len(taxids), batch_size), desc="Taxonomy batches"):
        batch = taxids[i:i + batch_size]
        params = {"db": "taxonomy", "id": ",".join(map(str, batch)), "retmode": "xml"}
        try:
            r = requests.post(ENTREZ_URL, data=params)
            r.raise_for_status()
            root = ET.fromstring(r.text)
            for taxon in root.findall("Taxon"):
                tid = taxon.find("TaxId").text
                family, order = "N/A", "N/A"
                lineage_ex = taxon.find("LineageEx")
                if lineage_ex is not None:
                    for entry in lineage_ex.findall("Taxon"):
                        rank = entry.find("Rank").text
                        name = entry.find("ScientificName").text
                        if rank == "family":
                            family = name
                        elif rank == "order":
                            order = name
                results[tid] = {"family": family, "order": order}
            time.sleep(0.35)  # Entrez allows 3 requests/s without an API key
        except Exception as e:
            print(f"Warning: taxonomy batch {i // batch_size + 1} failed: {e}")
    return results


def cmd_fetch(args):
    """All genome assemblies under TAXON_ID, with family and order, to OUTPUT_RAW."""
    print(f"Fetching assembly metadata for taxon {TAXON_ID}...")
    url = f"{API_BASE}/genome/taxon/{TAXON_ID}/dataset_report"
    params = {"page_size": 1000}
    all_genomes = []
    next_token = None
    page_num = 1

    try:
        while True:
            if next_token:
                params["page_token"] = next_token
            print(f"Page {page_num}...")
            response = requests.get(url, params=params)
            response.raise_for_status()
            data_json = response.json()
            reports = data_json.get('reports', [])
            if not reports and page_num == 1:
                sys.exit(f"No assemblies found for taxon {TAXON_ID}.")

            for record in reports:
                asm_info = record.get('assembly_info', {})
                org = record.get('organism', {})
                acc = record.get('accession') or record.get('current_accession')
                if not acc:
                    continue
                all_genomes.append({
                    "accession": acc,
                    "taxid": str(org.get('tax_id')),
                    "organism": org.get('organism_name'),
                    "assembly_level": asm_info.get('assembly_level', 'N/A'),
                    "submission_date": asm_info.get('release_date', '1900-01-01'),
                    "family": "N/A", "order": "N/A"
                })
            next_token = data_json.get('next_page_token')
            if not next_token:
                break
            page_num += 1
    except requests.RequestException as e:
        sys.exit(f"Error: {e}")

    unique_taxids = list(set(g['taxid'] for g in all_genomes if g['taxid']))
    taxonomy_map = get_lineage_batch(unique_taxids)
    for genome in all_genomes:
        tid = genome['taxid']
        if tid in taxonomy_map:
            genome['family'] = taxonomy_map[tid]['family']
            genome['order'] = taxonomy_map[tid]['order']

    pd.DataFrame(all_genomes).to_csv(OUTPUT_RAW, index=False)
    print(f"Saved {len(all_genomes)} assemblies to {OUTPUT_RAW}")


def cmd_filter(args):
    """Best `-n` assemblies per family (assembly level, then newest) to OUTPUT_FILTERED and ACCESSION_LIST."""
    if not Path(OUTPUT_RAW).exists():
        sys.exit(f"{OUTPUT_RAW} not found. Run 'fetch' first.")
    print(f"Keeping {args.number} assemblies per family...")
    df = pd.read_csv(OUTPUT_RAW)
    df = df[df['family'] != 'N/A']

    quality_map = {'Complete Genome': 3, 'Chromosome': 2, 'Scaffold': 1, 'Contig': 0}
    df['quality_score'] = df['assembly_level'].map(quality_map).fillna(0)
    df['submission_date'] = pd.to_datetime(df['submission_date'])

    df_sorted = df.sort_values(by=['quality_score', 'submission_date'], ascending=[False, False])
    selected_df = df_sorted.groupby('family').head(args.number)

    print(selected_df['family'].value_counts().head(5))
    selected_df.to_excel(OUTPUT_FILTERED, index=False)

    with open(ACCESSION_LIST, "w") as f:
        for acc in selected_df['accession']:
            f.write(f"{acc}\n")
    print(f"Selected {len(selected_df)} assemblies: {OUTPUT_FILTERED}, {ACCESSION_LIST}")


def cmd_download(args):
    """Genome FASTA, GFF and protein FASTA for ACCESSION_LIST, as one NCBI Datasets ZIP_FILE."""
    if not Path(ACCESSION_LIST).exists():
        sys.exit(f"{ACCESSION_LIST} not found. Run 'filter' first.")
    with open(ACCESSION_LIST) as f:
        accessions = [line.strip() for line in f if line.strip()]

    url = f"{API_BASE}/genome/download"
    payload = {
        "accessions": accessions,
        "include_annotation_type": ["GENOME_FASTA", "GENOME_GFF", "PROT_FASTA"]
    }

    try:
        print(f"Downloading {len(accessions)} genomes...")
        with requests.post(url, json=payload, stream=True) as r:
            r.raise_for_status()
            with open(ZIP_FILE, 'wb') as f:
                total_size = int(r.headers.get('content-length', 0))
                with tqdm(total=total_size, unit='B', unit_scale=True, desc="Download") as pbar:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
        print(f"Saved {ZIP_FILE}")
    except requests.RequestException as e:
        sys.exit(f"Error: {e}")


def cmd_organize(args):
    """Unzip ZIP_FILE into FINAL_DATASET_DIR/high_quality (complete, chromosome) and low_quality (scaffold, contig).

    FINAL_DATASET_DIR is deleted and rebuilt. Genomes not in OUTPUT_FILTERED are skipped.
    """
    if not Path(ZIP_FILE).exists():
        sys.exit(f"{ZIP_FILE} not found. Run 'download' first.")
    if not Path(OUTPUT_FILTERED).exists():
        sys.exit(f"{OUTPUT_FILTERED} not found. Run 'filter' first.")

    df = pd.read_excel(OUTPUT_FILTERED)
    acc_to_level = dict(zip(df['accession'], df['assembly_level']))

    temp_dir = Path("temp_ncbi_unzipped")
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    print("Unzipping...")
    with zipfile.ZipFile(ZIP_FILE, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    final_dir = Path(FINAL_DATASET_DIR)
    if final_dir.exists():
        shutil.rmtree(final_dir)
    hq_dir = final_dir / "high_quality"
    lq_dir = final_dir / "low_quality"
    hq_dir.mkdir(parents=True)
    lq_dir.mkdir()

    count_moved = 0
    for genome_dir in tqdm(list((temp_dir / "ncbi_dataset" / "data").iterdir())):
        if not genome_dir.is_dir():
            continue
        level = acc_to_level.get(genome_dir.name)
        if not level:
            continue

        level_lower = str(level).lower()
        target_folder = hq_dir if 'complete' in level_lower or 'chromosome' in level_lower else lq_dir
        for file_path in genome_dir.glob("*"):
            if file_path.is_file():
                shutil.copy(file_path, target_folder / file_path.name)
        count_moved += 1

    shutil.rmtree(temp_dir)
    print(f"Organised {count_moved} genomes into {FINAL_DATASET_DIR}/{{high_quality,low_quality}}")


def main():
    parser = argparse.ArgumentParser(description="Fetch, select and download cyanobacterial genomes from NCBI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("fetch", help=f"Assembly metadata and taxonomy -> {OUTPUT_RAW}")
    p_filter = subparsers.add_parser("filter", help=f"Best assemblies per family -> {OUTPUT_FILTERED}, {ACCESSION_LIST}")
    p_filter.add_argument("-n", "--number", type=int, default=10, help="Assemblies per family (default 10)")
    subparsers.add_parser("download", help=f"Download the selected genomes -> {ZIP_FILE}")
    subparsers.add_parser("organize", help=f"Unzip into {FINAL_DATASET_DIR}/ by assembly quality")

    args = parser.parse_args()
    {"fetch": cmd_fetch, "filter": cmd_filter, "download": cmd_download, "organize": cmd_organize}[args.command](args)


if __name__ == "__main__":
    main()
