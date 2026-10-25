"""Paths for the BGC annotation-gain analysis. Same GFF3 inputs as kegg_analyses/config.py."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

GENOME_LIST = ROOT / "data/raw/summaries/goodQ_summary_PBC.csv"  # `File` = <accession>_<asm>_genomic; also cross-check

ENHANCED_GFF_DIR = Path("/media/1GBDrive/cyano-annotator/Dataset_Good_Quality")
ENHANCED_GFF_GLOB = "{acc}_*/Reannota_New/enhanced/*_enhanced.gff3"  # one folder per genome
BAKTA_GFF_DIR = Path("/media/1GBDrive/cyano-annotator/benchmarks/cyanobacterial general/good_quality_reports/Bakta_gff")

OUT_DIR = HERE / "output"  # all working tables
EXPORT_DIR = ROOT / "data/raw/bgc_analyses"  # only the CSVs the manuscript uses

HYPOTHETICAL = "hypothetical protein"
# Same rule as kegg_analyses: these products replace "hypothetical protein" without naming a function (MobiDB-lite)
UNINFORMATIVE = {HYPOTHETICAL, "consensus disorder prediction"}
