"""Add Prokka annotation metrics to the per-genome summary CSV (edits it in place; also writes <csv>_enriched.xlsx)."""
import sys
import os
import logging
import re
from Bio import SeqIO
import pandas as pd

logger = logging.getLogger("ProkkaBatchUpdater")

def count_gff_features(gff_path):
    """Dbxref, hypothetical, pseudogene and BGC counts from a Prokka GFF file."""
    counts = {
        "Genes": 0, "GO": 0, "Go_richness": 0, "InterPro": 0, "PFAM": 0,
        "KEGG": 0, "COG": 0, "Pseudogene_candidates": 0, "BGCs": 0,
        "Hypothetical": 0, "EC": 0
    }

    if not gff_path or not os.path.exists(gff_path):
        return counts

    in_bgc = False
    gene_number = 0

    try:
        with open(gff_path, "r") as f:
            for line in f:
                if line.startswith("## Predicted BGC features"):
                    in_bgc = True
                    continue
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                if in_bgc:
                    if "biosynthetic-gene-cluster" in parts[2]:
                        counts["BGCs"] += 1
                    continue

                attributes = parts[8]
                counts["GO"] += attributes.count("GO:")
                counts["COG"] += attributes.count("COG:")
                counts["InterPro"] += attributes.count("InterPro:")
                counts["PFAM"] += attributes.count("PFAM:")
                counts["KEGG"] += attributes.count("KEGG")
                counts["EC"] += attributes.count("EC:") + attributes.count("eC_number=")

                if "product=hypothetical" in attributes:
                    counts["Hypothetical"] += 1

                if "Pseudogene candidate" in attributes or "pseudogene" in attributes:
                    counts["Pseudogene_candidates"] += 1

                if parts[2] in ['CDS', 'gene']:
                    gene_number += 1

        counts["Go_richness"] = (counts["GO"] / gene_number) if gene_number > 0 else 0

    except Exception as e:
        logger.error(f"Error parsing GFF {gff_path}: {e}")

    return counts

def get_gbk_metrics(gbk_file):
    """Described and undescribed (hypothetical) coding space, as % of bp, from a Prokka GenBank file."""
    if not gbk_file or not os.path.exists(gbk_file):
        return {'desc_pct': 0.0, 'undesc_pct': 0.0}

    total_described_bp = 0
    total_undescribed_bp = 0
    rna_types = ['tRNA', 'rRNA', 'ncRNA', 'tmRNA', 'mRNA', 'misc_RNA']

    try:
        for record in SeqIO.parse(gbk_file, "genbank"):
            rec_described = set()
            rec_undescribed = set()

            for feature in record.features:
                if feature.type in ['gene', 'source', 'assembly_gap']: continue

                current_indices = set()
                for part in feature.location.parts:
                    current_indices.update(range(int(part.start), int(part.end)))

                if feature.type == 'CDS':
                    product = feature.qualifiers.get('product', [''])[0].lower()
                    if 'hypothetical protein' in product:
                        rec_undescribed.update(current_indices)
                    else:
                        rec_described.update(current_indices)
                elif feature.type in rna_types:
                    rec_described.update(current_indices)

            total_described_bp += len(rec_described)
            total_undescribed_bp += len(rec_undescribed - rec_described)

    except Exception as e:
        logger.error(f"Error processing GBK '{gbk_file}': {e}")
        return {'desc_pct': 0.0, 'undesc_pct': 0.0}

    total_coding_space = total_described_bp + total_undescribed_bp
    if total_coding_space == 0:
        return {'desc_pct': 0.0, 'undesc_pct': 0.0}

    desc_pct = (total_described_bp / total_coding_space) * 100
    undesc_pct = (total_undescribed_bp / total_coding_space) * 100

    return {'desc_pct': desc_pct, 'undesc_pct': undesc_pct}

def extract_accession_id(filename_string):
    """GCA_/GCF_ accession from a CSV File value."""
    match = re.search(r'(GC[AF]_\d+)', str(filename_string))
    if match:
        return match.group(1)
    return None

def find_prokka_file(search_id, folder, extension):
    if not folder or not os.path.exists(folder):
        return None
    for f in os.listdir(folder):
        if f.endswith(extension) and search_id in f:
            return os.path.join(folder, f)
    return None

def update_all_genomes(csv_path, prokka_gff_dir, prokka_gbk_dir):
    logger.info(f"Loading summary CSV: {csv_path}")
    if not os.path.exists(csv_path):
        logger.error("CSV file not found.")
        sys.exit(1)

    df = pd.read_csv(csv_path)

    if 'File' not in df.columns:
        logger.error("CSV must have a 'File' column containing the Genome IDs.")
        sys.exit(1)

    logger.info(f"Found {len(df)} genomes to process.")

    for index, row in df.iterrows():
        genome_id_full = row['File']
        search_id = extract_accession_id(genome_id_full)
        if not search_id:
            search_id = str(genome_id_full).split('_')[0]

        prokka_gff = find_prokka_file(search_id, prokka_gff_dir, ".gff")
        if not prokka_gff: prokka_gff = find_prokka_file(search_id, prokka_gff_dir, ".gff3")
        prokka_gbk = find_prokka_file(search_id, prokka_gbk_dir, ".gbk")

        if not prokka_gff:
            logger.warning(f"[{search_id}] No Prokka GFF found. Skipping metrics.")
            continue

        p_gff_counts = count_gff_features(prokka_gff)
        p_gbk_metrics = get_gbk_metrics(prokka_gbk) if prokka_gbk else {'desc_pct': 0, 'undesc_pct': 0}

        df.at[index, "Prokka_GO_entries"] = p_gff_counts["GO"]
        df.at[index, "Prokka_GO_richness"] = p_gff_counts["Go_richness"]
        df.at[index, "Prokka_InterPro_entries"] = p_gff_counts["InterPro"]
        df.at[index, "Prokka_PFAM_entries"] = p_gff_counts["PFAM"]
        df.at[index, "Prokka_COG_entries"] = p_gff_counts["COG"]
        df.at[index, "Prokka_KEGG_entries"] = p_gff_counts["KEGG"]
        df.at[index, "Prokka_Pseudogene_candidates"] = p_gff_counts["Pseudogene_candidates"]
        df.at[index, "Prokka_BGCs"] = p_gff_counts["BGCs"]
        df.at[index, "Prokka_EC_entries"] = p_gff_counts["EC"]
        df.at[index, "Prokka_Hypotheticals"] = p_gff_counts["Hypothetical"]
        df.at[index, "Prokka_Desc_PCT"] = f"{p_gbk_metrics['desc_pct']:.2f}%"
        df.at[index, "Prokka_Undesc_PCT"] = f"{p_gbk_metrics['undesc_pct']:.2f}%"

    try:
        df.to_csv(csv_path, index=False)
        logger.info(f"Success! Updated CSV saved to: {csv_path}")
    except Exception as e:
        logger.error(f"Failed to save CSV: {e}")

    for col in ['Prokka_Desc_PCT', 'Prokka_Undesc_PCT']:
        if col in df.columns:
            df[col] = df[col].str.replace('%', '').astype(float)

    excel_out = csv_path.replace('.csv', '_enriched.xlsx')
    df.to_excel(excel_out, index=False)
    logger.info(f"Excel version for R saved at: {excel_out}")

