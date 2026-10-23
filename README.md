# manuscript-code

Analysis code for the manuscript on cyanobacterial genome annotation enrichment. The evaluated pipeline (Cyan-ote) runs Bakta, then InterProScan and eggNOG-mapper, and merges the results with Reannota. It is compared with Bakta alone and with Prokka.

| Folder | Analysis |
|---|---|
| [genome_analyzer/](genome_analyzer/) | Per-genome statistics and plots: hypothetical-protein reduction, Dbxref gains, Prokka metrics, antiSMASH vs GECCO clusters |
| [kegg_analyses/](kegg_analyses/) | KO coverage and KEGG module reconstruction, Bakta vs integrated pipeline, 40 genomes |
| [bgc_analyses/](bgc_analyses/) | Annotation gain for CDS inside antiSMASH clusters |

## Setup

Python 3.13.

```
pip install -r requirements.txt
```

## Inputs

Not included in this repository. Input paths are set in `kegg_analyses/config.py` and `bgc_analyses/config.py`:

- Bakta GFF3s (`BAKTA_GFF_DIR`) and Reannota enhanced GFF3s (`ENHANCED_GFF_DIR`), one per genome, named `<accession>_*.gff3`.
- Under `data/` at the repo root: `raw/summaries/Genomes.xlsx` (genome manifest), `processed/GOODQ_PBC_Summary.csv` (per-genome pipeline summary) and `processed/BGC_statistics_summary.xlsx` (cluster counts, cross-check only).

The KEGG and InterPro reference snapshot used for the manuscript is included in `kegg_analyses/reference/`.

## Run

All commands run from the repo root. See each folder's README.

```
python -m pytest kegg_analyses/tests
python -m pytest bgc_analyses/tests
```

Run the two test suites separately, because both folders have a module named `config`.
