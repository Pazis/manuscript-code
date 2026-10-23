"""Paths, module panel and expectations for the KEGG analysis."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

GENOMES_XLSX = ROOT / "data/raw/summaries/Genomes.xlsx"
MANIFEST_SHEET = "KO Analysis Genomes"
PIPELINE_SUMMARY = ROOT / "data/processed/GOODQ_PBC_Summary.csv"

ENHANCED_GFF_DIR = Path("/media/1GBDrive/cyano-annotator/GoodQ_GFF")
BAKTA_GFF_DIR = Path("/media/1GBDrive/cyano-annotator/benchmarks/cyanobacterial general/good_quality_reports/Bakta_gff")

KEGG_DIR = HERE / "reference/kegg"
PFAM2IPR = HERE / "reference/interpro/pfam2interpro.tsv"  # eggNOG Pfam names -> InterPro, for IPS confirmation
OUT_DIR = HERE / "output"  # all working tables
EXPORT_DIR = ROOT / "data/raw/kegg_analyses"  # only the CSVs the manuscript uses (supplementary, data availability, plots)
TABLES_DIR = ROOT / "tables"  # manuscript Tables 1 and 2

# Manifest genomes without good-quality outputs -> replacements
SUBSTITUTIONS = {
    "GCF_000252425.1": ("GCA_965198675.1", "Prochloron sp.",
                        "Replaces Prochloron didemni P2-Fiji (scaffold, no good-quality outputs); same family"),
    "GCF_000332315.1": ("GCF_000018105.1", "Acaryochloris marina MBIC11017",
                        "Replaces Prochlorothrix hollandica PCC 9006 (no Prochlorotrichaceae genome with full outputs); Chl d antenna"),
}

PANEL = ["M00161", "M00163", "M00165", "M00004", "M00854", "M00175", "M00531", "M00611"]
COMPOSITE = {"M00611"}
BENCHMARK_ORGS = {"GCF_050408035.1": "syz", "GCF_000025125.1": "cyu"}  # Synechocystis 6803, UCYN-A

# Expected completeness for UCYN-A (Candidatus Atelocyanobacterium thalassa ALOHA); False = must not be complete
EXPECTATIONS = {
    "GCF_000025125.1": {"M00161": False, "M00165": False, "M00611": False, "M00163": True, "M00175": True},
}

# T2 tier: KEGG Metabolism maps, minus the global/overview maps (map011xx, map012xx)
def is_metabolic_map(map_id: str) -> bool:
    n = int(map_id.removeprefix("map"))
    return 10 <= n < 1100

HYPOTHETICAL = "hypothetical protein"
