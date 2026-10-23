# BGC analyses

**Question:** of the CDS inside antiSMASH clusters, how many are better annotated by the integrated pipeline (Bakta → InterProScan + eggNOG-mapper → Reannota, "Cyan-ote") than by Bakta alone? Measured as hypothetical proteins recovered and Dbxrefs gained. GECCO clusters and CDS outside clusters are not analysed.

**Inputs** (same GFF3s as `kegg_analyses/`; paths in `config.py`):
- Genomes: every `File` accession in `data/processed/GOODQ_PBC_Summary.csv` with exactly one GFF3 in each directory. Others are skipped and listed.
- Baseline (state A): Bakta GFF3s.
- Integrated (state D): Reannota enhanced GFF3s. These also hold the clusters, as `biosynthetic-gene-cluster` features with source `antiSMASH[:version]`.
- Cross-check: `Enhanced_BGCs` in `GOODQ_PBC_Summary.csv` and `AntiSMASH Count` in `data/processed/BGC_statistics_summary.xlsx`.

**Run** (from the repo root):

```
python bgc_analyses/run_bgc_eval.py all      # or build | export
python -m pytest bgc_analyses/tests
```

## Definitions

- **Cluster:** an antiSMASH protocluster as written by Reannota, not an antiSMASH region. Overlapping protoclusters (e.g. NRPS-like inside NRPS) are separate features.
- **In a cluster:** the Bakta CDS lies entirely within a cluster on the same contig. The CDS copies Reannota writes under each cluster are ignored. A CDS in overlapping protoclusters is counted once; `bgcs` lists all of them.
- **Hypothetical / recovered:** as in `kegg_analyses`. The Bakta product is "hypothetical protein", and the enhanced product is neither that nor "consensus disorder prediction" (MobiDB-lite). Counts are per CDS.
- **Dbxrefs:** distinct `KEGG`, `GO`, `COG` and `EC` entries in the CDS's `Dbxref` attribute; all other databases are ignored. COG category letters (eggNOG `COG_cat:<letter>`, Bakta `COG:<letter>`) are dropped, so only COG IDs (`COG:COG####`) count. `EC:EC:` is normalised to `EC:`. Added = in enhanced but not in Bakta; lost = the reverse.
- **Improved:** recovered, or at least one Dbxref added.

## Outputs

Working tables go to `output/` (TSV). `export` copies all except the cross-check to `data/raw/bgc_analyses/` (CSV).

| File | Rows |
|---|---|
| `summary` | One row: totals, pooled %, and median [IQR] of per-genome % |
| `per_genome` | One per genome (zeros if no cluster): cluster, CDS, hypothetical, recovered, gained-Dbxref and improved counts; Dbxref totals |
| `dbxref_by_db` | Genome × database (KEGG, GO, COG, EC): Dbxref entries in Bakta, in enhanced, and added |
| `bgc_cds` | One per CDS inside a cluster: cluster IDs and types, products before and after, gains |
| `clusters` | One per antiSMASH cluster: ID, contig, coordinates, product |
| `crosscheck_cluster_counts` | Parsed clusters vs `Enhanced_BGCs` and `AntiSMASH Count`. A mismatch with the workbook may reflect protoclusters vs regions. |

## Caveat

eggNOG adds GO ancestor terms, so GO dominates "Dbxrefs added". In a trial on the Reannota example genome, 116 of the 280 cluster CDS that gained a Dbxref gained only GO terms. Report gains per database using `dbxref_by_db` and `added_dbs`.
