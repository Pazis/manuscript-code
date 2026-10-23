# KEGG analyses

**Question:** what metabolic annotation does the integrated pipeline (Bakta → InterProScan + eggNOG-mapper → Reannota) add to the Bakta baseline in 40 cyanobacteria? Measured as KO coverage and KEGG module reconstruction.

**Inputs** (paths in `config.py`):
- Manifest: `data/raw/summaries/Genomes.xlsx`, sheet `KO Analysis Genomes`, with the 2 substitutions in `config.SUBSTITUTIONS`.
- Baseline (state A): Bakta GFF3s.
- Integrated (state D): Reannota enhanced GFF3s.
- Reference snapshot: `reference/kegg/` (KEGG Module release 2026/10/09) and `reference/interpro/pfam2interpro.tsv` (InterPro 110.0).
- Cross-check: `data/processed/GOODQ_PBC_Summary.csv`.

**Run** (from the repo root):

```
python kegg_analyses/run_kegg_eval.py all      # or validate | build | reconstruct | trace | review | export | table1
python kegg_analyses/run_kegg_eval.py table2   # after curating data/raw/kegg_analyses/review_candidates.csv
python -m pytest kegg_analyses/tests
```

`fetch_kegg_reference.py` rebuilds `reference/` from the current KEGG and InterPro releases. Do not run it to reproduce the manuscript.

**Outputs:**
- `output/`: all working tables, described in [output/README.md](output/README.md).
- `data/raw/kegg_analyses/`: the CSVs used in the manuscript, written by `export`.
- `tables/kegg_table1.csv`, `tables/kegg_table2.csv`: Tables 1 and 2.

## Method notes

The inputs are GFF3s, not raw tool outputs. This determines the following:

- **Provenance:** a KO on the enhanced CDS that Bakta did not assign is credited to eggNOG-mapper, because InterProScan assigns no KOs. `validate` fails if any Bakta CDS or Bakta KO is missing from the enhanced GFF3, since gains are computed as D − A.
- **Bakta and eggNOG agreement** cannot be observed, because the GFF3s don't record a KO that both tools gave. `provenance` is therefore `Bakta` or `eggNOG`.
- **Conflicts:** Reannota adds KOs only to proteins without a Bakta KO, at most one per protein. The `conflict` column is 0 everywhere.
- **eggNOG fields:** e-value, score, `max_annot_lvl` and multi-KO rows are not in the GFF3s, so there is no ambiguity flag and no multi-KO review criterion.
- **IPS confirmation:** eggNOG Pfam names (`Note: Eggnog PFAM comment:`) are mapped to InterPro entries and checked for overlap with the protein's `InterPro:` Dbxrefs.
- **Hypothetical recovery:** the enhanced product is neither "hypothetical protein" nor "consensus disorder prediction" (MobiDB-lite).
- **antiSMASH features:** the enhanced GFF3s repeat BGC genes as extra `antiSMASH` CDS features. These are skipped.
- **Tool versions** are not recorded in the GFF3s. `provenance.tsv` records each file's md5 and header instead.
- **T2 tier:** KOs in maps map00010–map01099, which excludes the global and overview maps (map011xx, map012xx).
- **M00611** is a composite module that KEGG does not link to organisms, so it never counts as a disagreement with KEGG.
- **Low Bakta baseline:** no panel module is complete in the baseline for any genome. For Synechocystis 6803, Bakta's own `.tsv` (v1.11.4, db v6.0) has the same 725 KOs as its GFF3. Bakta often names a protein correctly without assigning a KO.
- **Requires review:** disagreement with KEGG's organism annotation is checked only for Synechocystis (`syz`) and UCYN-A (`cyu`), each against its own KEGG code.
- **Review candidates:** at most 3 per genome and 2 per KO, so that repeats of one case don't dominate the list.
