# KEGG analysis outputs

Written by `kegg_analyses/run_kegg_eval.py`; do not edit by hand. Method notes are in [../README.md](../README.md).

All files are tab-separated. `genome` is always the assembly accession. `protein` is the CDS ID shared by the Bakta and enhanced GFF3s (e.g. `OMICDL_00005`).

**Annotation states**
- **A**: Bakta baseline.
- **B**: Bakta + InterProScan. Equal to A for KO metrics, because InterProScan assigns no KOs.
- **C**: Bakta ∪ eggNOG-mapper.
- **D**: full integrated pipeline (Reannota enhanced GFF3). Equal to C for KO metrics.

**Analysis sets**
- `all`: every Bakta CDS.
- `hypothetical`: CDS whose Bakta product is "hypothetical protein".

**KO tiers** (fixed before any results were viewed)
- **T1**: KOs that appear in the 8 panel module definitions (M00161, M00163, M00165, M00004, M00854, M00175, M00531, M00611).
- **T2**: KOs linked to KEGG Metabolism maps map00010–map01099. The global/overview maps map011xx and map012xx are excluded.

---

## Inputs and validation (`validate`)

### `genomes.tsv` — 40 rows, one per genome
The 40-genome manifest: the `KO Analysis Genomes` sheet of `data/raw/summaries/Genomes.xlsx` with 2 substitutions applied.

| Column | Meaning |
|---|---|
| `group`, `accession`, `organism`, `family`, `assembly_level`, `role` | From the sheet. |
| `substitution_note` | Non-empty for the 2 replaced genomes: GCA_965198675.1 replaces *Prochloron didemni*, and GCF_000018105.1 replaces *Prochlorothrix hollandica*. |
| `kegg_org` | KEGG organism code for the benchmark genomes (`syz` = Synechocystis 6803, `cyu` = UCYN-A); otherwise empty. |
| `bakta_gff`, `enhanced_gff` | Absolute paths of the input GFF3s. |

### `provenance.tsv` — 120 rows, 3 per genome
Fingerprints of the input files. Tool versions are not recorded in the GFF3s.

| Column | Meaning |
|---|---|
| `file` | One of: `bakta` (Bakta GFF3), `enhanced` (enhanced GFF3), or `enhanced_only_cds` (count of CDS IDs found only in the enhanced GFF3, after antiSMASH features are skipped). |
| `path` | Input file path. |
| `n_cds` | CDS count, or the enhanced-only count for `enhanced_only_cds` rows. |
| `md5` | Checksum of the input file. |
| `header` | The GFF3's `#` header lines. |

---

## Protein-level evidence (`build`)

### `proteins.tsv` — 199,685 rows, one per Bakta CDS
Per-protein annotation from both GFF3s. List fields are comma-separated.

| Column | Meaning |
|---|---|
| `bakta_product`, `bakta_gene` | Bakta product and gene name. |
| `enhanced_product`, `enhanced_name` | Product and name in the enhanced GFF3. |
| `is_hypothetical` | The Bakta product is "hypothetical protein". |
| `bakta_kos`, `enhanced_kos` | KOs in state A and in state D. |
| `new_kos` | `enhanced_kos − bakta_kos`, i.e. the KOs added by eggNOG-mapper. |
| `dropped_kos` | `bakta_kos − enhanced_kos`. Always empty, because validation fails otherwise. |
| `bakta_ec`, `enhanced_ec` | EC numbers in each GFF3. |
| `interpro` | InterPro IDs on the enhanced CDS (from InterProScan). |
| `eggnog_pfams` | Pfam names reported by eggNOG-mapper (`Eggnog PFAM comment` in the Note). |
| `has_ips` | Any InterProScan evidence: InterPro IDs, or an "InterPro detected" note. |
| `has_eggnog` | Any eggNOG-mapper evidence: Pfam names, or a new KO. |
| `ips_confirmed` | eggNOG's Pfams, mapped to InterPro entries, overlap the protein's InterProScan entries. |
| `named_after` | The enhanced product is neither "hypothetical protein" nor "consensus disorder prediction". |

### `protein_ko.tsv` — 84,506 rows, one per (genome, protein, KO) in state D
The protein–KO associations, with provenance.

| Column | Meaning |
|---|---|
| `in_bakta` | Bakta assigned this KO, so the pair is in state A. |
| `in_eggnog` | The KO was added by eggNOG-mapper. This is always `not in_bakta`, because agreement between the two tools cannot be observed in the GFF3s. |
| `provenance` | `Bakta` or `eggNOG`. |
| `ips_confirmed`, `is_hypothetical` | Copied from `proteins.tsv`. |
| `conflict` | eggNOG added a KO to a protein that already had a different Bakta KO. Always False with this pipeline. |

---

## Coverage metrics (`build`)

### `coverage_per_state.tsv` — 640 rows (40 genomes × 2 sets × 2 tiers × 4 states)

| Column | Meaning |
|---|---|
| `N` | Number of proteins in the set. |
| `M` | Proteins with ≥1 qualifying KO in this state. |
| `pct` | `100·M/N`. |
| `pairs` | Qualifying protein–KO pairs. |
| `distinct_kos` | Distinct qualifying KO IDs. |
| `conflicts` | Number of conflict pairs. |

### `coverage_per_genome.tsv` — 160 rows (40 genomes × 2 sets × 2 tiers)
Baseline (A) and integrated (D) coverage side by side.

| Column | Meaning |
|---|---|
| `M_baseline`, `M_integrated` | Proteins with ≥1 qualifying KO in A and in D. |
| `pct_baseline`, `pct_integrated` | The same, as a percentage of `N`. |
| `gain_pp` | Gain in percentage points. |
| `new_pairs` | Protein–KO pairs in D but not in A. |
| `new_distinct_kos` | KO IDs in D but not in A. |
| `conflicts` | As above. |

Exported for Figure 8 and Figure S1 (the paired A-vs-D plot).

### `hypothetical_recovery.tsv` — 40 rows, one per genome
Counts of Bakta-hypothetical proteins.

| Column | Meaning |
|---|---|
| `n_hypothetical` | Total Bakta-hypothetical proteins. |
| `has_ips` | With any InterProScan evidence. |
| `has_eggnog` | With any eggNOG-mapper evidence. |
| `any_signal` | With either. |
| `named_after` | Given an informative product in the enhanced GFF3. |
| `new_T1_ko` | With a new T1 KO. |
| `new_T2_ko` | With a new T2 KO. |
| `new_ec` | With a new EC number. |

### `table1_global.tsv` — 44 rows
Draft of Table 1: summaries across the 40 genomes.

| Column | Meaning |
|---|---|
| `metric` | A metric name from `coverage_per_genome.tsv`, or `hypothetical_*` for the columns of `hypothetical_recovery.tsv`. |
| `set`, `tier` | Analysis set and KO tier (`tier` is empty for the `hypothetical_*` rows). |
| `n_genomes` | Number of genomes summarised. |
| `median`, `q1`, `q3` | Median and interquartile range across genomes. |
| `total` | Sum over genomes, for count metrics only. |
| `wilcoxon_p` | Paired Wilcoxon signed-rank p-value comparing `pct_integrated` with `pct_baseline`. Given on the `gain_pp` rows only. |

### `crosscheck_pipeline_summary.tsv` — 40 rows
A sanity check against `data/processed/GOODQ_PBC_Summary.csv`.

| Column | Meaning |
|---|---|
| `distinct_A`, `distinct_D` | Distinct KOs per genome in A and in D. |
| `Initial_KEGG_entries`, `Enhanced_KEGG_entries` | The pipeline's Dbxref occurrence counts. |
| `ratio_A`, `ratio_D` | Occurrences ÷ distinct KOs. |

The ratios are expected to be ≥1 and of the same order of magnitude, not equal to 1, because the pipeline counts occurrences rather than distinct KOs. Observed range: 1.05–1.87.

---

## Module reconstruction (`reconstruct`)

### `modules_per_genome.tsv` — 1,280 rows (40 genomes × 4 states × 8 modules)
Module completeness, following the KEGG Mapper block convention.

| Column | Meaning |
|---|---|
| `steps_total` | Number of top-level blocks in the module definition. |
| `steps_present` | Blocks satisfied. |
| `blocks_missing` | Blocks not satisfied. |
| `complete` | True when no block is missing. |
| `status` | One of `complete`, `1 block missing`, `2 blocks missing`, `incomplete`. |

M00611 is a composite of M00161, M00163 and M00165. Do not count it as an independent gain.

### `kegg_reference_check.tsv` — 16 rows (2 organisms × 8 modules)
Checks the module evaluator against KEGG's own annotations of `syz` and `cyu`.

| Column | Meaning |
|---|---|
| `status` | Our evaluator's result on KEGG's KO set for that organism. |
| `kegg_lists_module` | Whether KEGG links the module to the organism. |
| `agree` | Our result matches KEGG's. M00611 always counts as agreeing, because KEGG does not link composite modules to organisms. |

All 16 rows agree.

---

## Traceability (`trace`)

### `supp_genome_module.tsv` — 320 rows (40 genomes × 8 modules)
Table S3: module status in A and D, with the eggNOG KOs behind each change.

| Column | Meaning |
|---|---|
| `baseline_status`, `integrated_status` | Status in A and D, as in `modules_per_genome.tsv`. |
| `baseline_steps`, `integrated_steps` | Blocks present / total. |
| `change_category` | `complete in both`, `newly complete`, `improved` (more blocks in D) or `unchanged`. |
| `contributing_kos`, `contributing_proteins` | New eggNOG KOs that satisfy blocks satisfied only in D (`protein:KO`). For M00611, those of its child modules. |
| `n_ips_confirmed`, `n_hypothetical_in_bakta` | Contributing pairs that are IPS-confirmed, and that were Bakta-hypothetical. |
| `requires_review`, `review_reasons` | Flagged when a change rests only on eggNOG without IPS confirmation, contradicts `config.EXPECTATIONS`, or disagrees with KEGG's organism annotation (`syz`, `cyu` only). |
| `interpretation` | Empty; filled in by hand. |

### `module_contributions.tsv` — one row per (genome, module, contributing protein–KO pair)
The long form of `contributing_proteins`. Columns: `ko`, `protein`, `bakta_product`, `is_hypothetical`, `ips_confirmed`, `change_category`, `requires_review`. `review` uses it to find which new KOs affect a module.

### `nitrogen_table.tsv` — 40 rows
Table S4. One column per nitrogen-metabolism gene (`traceability.NITROGEN_KOS`), listing `protein(provenance)` hits; Bakta and enhanced *nif/vnf/anf* gene names; M00175 and M00531 status in A and D, with their new contributors.

### `benchmark_syz_cyu.tsv` — 16 rows (2 genomes × 8 modules)
Table S5. KO counts and Jaccard similarity to KEGG's own KO set for `syz` and `cyu`, in A and D; module status on KEGG's KOs, in A and in D.

### `tca_new_assignments.tsv` — new eggNOG KOs in the TCA cycle map (map00020)
A qualitative list only, with no pass/fail. Columns: `genome`, `protein`, `ko`, `ko_name`, `bakta_product`, `is_hypothetical`, `ips_confirmed`.

## Review (`review`)

### `review_candidates_new.csv`
Written instead of `data/raw/kegg_analyses/review_candidates.csv` when that curated file already exists, which is never overwritten.
