# genome_analyzer

Per-genome annotation statistics and plots. Run from the repo root:

```
python -m genome_analyzer.main stats  -i <summary.csv> -g <gff_dir> -o plots
python -m genome_analyzer.main bgc    -d <enhanced_gff_dir> -o bgc_results
python -m genome_analyzer.main prokka --csv_file <summary.csv> --prokka_gff_dir <dir> --prokka_gbk_dir <dir>
```

Each command writes `<command>.log` to its output directory.

| Command | Input | Output |
|---|---|---|
| `stats` | Per-genome summary CSV (`File`, `Initial_*`, `Enhanced_*`, optional `Prokka_*` columns). `-g` adds contig count and genome size from the GFFs. | Figures 1–7 and `pipeline_stats_summary.xlsx` |
| `bgc` | Directory of enhanced GFF3s with antiSMASH and GECCO `biosynthetic-gene-cluster` features | `bgc_analysis_results.xlsx`, `bgc_statistical_summary.xlsx`, Figures 1–4 (`--skip-plots` skips them) |
| `prokka` | Summary CSV and Prokka GFF and GBK directories | Adds `Prokka_*` columns to the CSV **in place**, and writes `<csv>_enriched.xlsx` |

## Methods

- **Undescribed coding space (Prokka):** bp of "hypothetical protein" CDS as a percentage of bp in all CDS and RNA features, overlaps counted once. The `Initial_*` (Bakta) and `Enhanced_*` values come from the pipeline summary.
- **Initial vs Enhanced:** paired t-test if the per-genome differences pass Shapiro–Wilk (p > 0.05), otherwise Wilcoxon signed-rank. The same rule applies to antiSMASH vs GECCO cluster counts.
- **Genome size vs hypothetical %:** Pearson if both variables pass Shapiro–Wilk, otherwise Spearman.
- **BGC overlap:** clusters overlap when their coordinates intersect. Consensus = overlapping antiSMASH clusters / all antiSMASH clusters. Efficiency = 2 × overlap bp / (antiSMASH bp + GECCO bp). Consolidation ratio = antiSMASH clusters / GECCO clusters. Statistics use only genomes where both tools found overlapping clusters.
