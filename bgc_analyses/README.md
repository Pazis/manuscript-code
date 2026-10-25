# BGC analyses

**Question:** Of the CDS inside antiSMASH-defined biosynthetic gene clusters, how many are better annotated by the integrated pipeline (Bakta → InterProScan + eggNOG-mapper → Reannota, "Cyan-ote") than by the Bakta baseline? Measured as hypothetical proteins recovered and Dbxrefs gained. GECCO clusters and CDS outside clusters are not analysed.

**Inputs** (same GFF3s as `kegg_analyses/`, paths in `config.py`):
- Genome list: every accession in the `File` column of `data/raw/summaries/goodQ_summary_PBC.csv` that has exactly one GFF3 in each directory (others are skipped and listed).
- Baseline (state A): Bakta GFF3s in `/media/1GBDrive/cyano-annotator/benchmarks/cyanobacterial general/good_quality_reports/Bakta_gff`.
- Integrated (state D): Reannota enhanced GFF3s in `/media/1GBDrive/cyano-annotator/Dataset_Good_Quality/<genome>/Reannota_New/enhanced/`. These also hold the clusters (`biosynthetic-gene-cluster` features with source `antiSMASH[:version]`; `GECCO` ones are ignored).
- Cross-check: `Enhanced_BGCs` in the same genome list.

**How to run** (from the repo root):

```
python bgc_analyses/run_bgc_eval.py all      # or build | export
python -m pytest bgc_analyses/tests
```

## Definitions

- **Cluster:** an antiSMASH **region** (`as:<contig>_region<n>`), which is what Reannota_New writes. The old GoodQ_GFF run wrote overlapping protoclusters instead.
- **In a cluster:** the Bakta CDS lies entirely inside a cluster's coordinates on the same contig. The CDS copies that Reannota writes under each cluster are ignored. A CDS in overlapping clusters is counted once, and `bgcs` lists all of them.
- **Hypothetical / recovered:** same rule as `kegg_analyses`. The Bakta product is "hypothetical protein", and the enhanced product is neither that nor "consensus disorder prediction" (MobiDB-lite). Counts are in CDS space.
- **Dbxrefs:** distinct `KEGG`, `GO`, `COG` and `EC` entries in the CDS's `Dbxref` attribute. Everything else is ignored, including InterPro, SO, RefSeq and UniRef. Two kinds of COG category letter are dropped: eggNOG's `COG_cat:<letter>` and Bakta's `COG:<letter>`. Only COG IDs (`COG:COG####`) count. `EC:EC:` is normalised to `EC:` so it isn't counted as a gain. Added = in enhanced but not in Bakta; lost = the reverse.
- **Improved:** recovered, or at least one Dbxref added.

## Outputs

`output/` (working tables, TSV) → `data/raw/bgc_analyses/` (CSV, written by `export`):

| File | Rows |
|---|---|
| `summary` | One row: totals across genomes, pooled %, and median [IQR] of per-genome % |
| `per_genome` | Per genome (zeros if no cluster): clusters, CDS, hypotheticals, recovered, gained-Dbxref and improved counts, Dbxref totals |
| `dbxref_by_db` | Genome × database (KEGG, GO, COG, EC): Dbxref entries in Bakta, enhanced, added |
| `bgc_cds` | Every CDS inside a cluster, with its cluster IDs and types, products before/after and gains |
| `clusters` | Every antiSMASH cluster parsed: ID, contig, coordinates, product |
| `crosscheck_cluster_counts` (output only) | Parsed clusters vs `Enhanced_BGCs` in the genome list |

## Caveats

- eggNOG adds GO ancestor terms, so GO dominates "Dbxrefs added". In a trial on the ReAnnota example genome, 116 of 280 cluster CDS that gained a Dbxref gained only GO terms. Use `dbxref_by_db` and `added_dbs` to report gains per database.

## Consolidation and boundary overlap (Table S8, Figures 6–7)

`bgc_consolidation.py` compares antiSMASH and GECCO cluster by cluster in the Reannota enhanced GFF3s
(`/media/1GBDrive/cyano-annotator/Dataset_Good_Quality/*/Reannota_New/enhanced/*_enhanced.gff3`, 164 genomes).

```
python bgc_analyses/bgc_consolidation.py                 # TSVs in output/ + sheets "S8 BGC consolidation" and "S8 summary"
python bgc_analyses/bgc_consolidation.py --no-workbook   # TSVs only
Rscript figures/fig06_bgc_consolidationratio/bgc_consolidation.r
Rscript figures/fig07_bgc_boundary_precision/boundary_precision.r
```

- **Unit:** one antiSMASH **region** vs one GECCO cluster. GFF3s from before the ReAnnota fix (protocluster IDs
  `as:<contig>_bgc<n>`) are refused, because overlapping protoclusters inflated the antiSMASH counts (2,342 protoclusters
  vs 1,673 regions in this cohort) and the overlap figures.
- **Overlap:** clusters overlap only on the same contig (≥ 1 bp). The previous table compared coordinates across contigs.
- **Overlap efficiency:** mean Jaccard index (shared bp / union bp) over overlapping region–cluster pairs, so it lies in
  0–1 (six genomes in the previous table were above 1).
- **Checks** printed on every run: efficiency and consensus in 0–1, coverage in 0–100 %, matched counts ≤ cluster counts,
  overlaps found from both sides, no overlapping antiSMASH regions, clusters only on annotated contigs, duplicate accessions.
- Column definitions are in the script's docstring. Outputs: `output/bgc_consolidation_per_genome.tsv` (= sheet S8),
  `output/bgc_consolidation_pairs.tsv` (every overlapping region–cluster pair), `output/bgc_consolidation_summary.tsv`
  (= sheet "S8 summary").
