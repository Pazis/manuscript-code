# genome_fetcher

Selects and downloads cyanobacterial genomes (NCBI taxon 1117) from NCBI Datasets. Run the steps in order, in the directory where the files should go:

```
python genome_fetcher/main.py fetch           # all assemblies + family/order -> cyanobacteria_raw.csv
python genome_fetcher/main.py filter -n 10    # best n per family -> cyanobacteria_selected.xlsx, accession_list.txt
python genome_fetcher/main.py download        # genome FASTA, GFF, protein FASTA -> genomes_dataset.zip
python genome_fetcher/main.py organize        # -> Cyano_Benchmark_Dataset/{high_quality,low_quality}/
```

- **Taxonomy:** family and order come from NCBI Taxonomy (Entrez). Assemblies with no family are dropped by `filter`.
- **Selection:** within each family, assemblies are ranked by assembly level (Complete Genome > Chromosome > Scaffold > Contig), then by release date, newest first.
- **Quality folders:** `high_quality` holds Complete Genome and Chromosome assemblies; `low_quality` holds Scaffold and Contig.
- `organize` deletes and rebuilds `Cyano_Benchmark_Dataset/`.

NCBI content changes over time, so rerunning `fetch` will not reproduce the original selection. The genomes used are listed in `data/raw/summaries/Genomes.xlsx`, sheet `Cyanobacterial Genomes`.
