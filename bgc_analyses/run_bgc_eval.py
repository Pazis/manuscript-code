"""BGC annotation gain: are CDS inside antiSMASH clusters better annotated after Reannota?

Run from the repo root:
  python bgc_analyses/run_bgc_eval.py {build,export,all}
Working tables go to bgc_analyses/output/; export copies the manuscript subset as CSV to data/raw/bgc_analyses/.
"""
import sys

import pandas as pd

import bgc
import config

OUT = config.OUT_DIR


def write(df: pd.DataFrame, name: str):
    df.to_csv(OUT / name, sep="\t", index=False)
    print(f"  wrote {name} ({len(df)} rows)")


def find_one(directory, accession, pattern="{acc}_*.gff3") -> str:
    hits = sorted(directory.glob(pattern.format(acc=accession)))
    return str(hits[0]) if len(hits) == 1 else ""


def genome_list() -> pd.DataFrame:
    df = pd.read_csv(config.GENOME_LIST)
    df["File"] = df.File.str.extract(r"^(GC[AF]_\d+\.\d+)", expand=False)  # folder name -> accession
    assert df.File.notna().all() and df.File.is_unique, "genome list has unparsable or duplicate accessions"
    return df


def build():
    genomes = sorted(genome_list().File)
    cds, clusters, dbs, skipped, analysed = [], [], [], [], []
    for acc in genomes:
        b, e = find_one(config.BAKTA_GFF_DIR, acc), find_one(config.ENHANCED_GFF_DIR, acc, config.ENHANCED_GFF_GLOB)
        if not b or not e:
            skipped.append(f"{acc} (bakta={bool(b)}, enhanced={bool(e)})")
            continue
        c, cl = bgc.cds_table(acc, b, e)
        analysed.append(acc)
        dbs.append(bgc.per_db(c))
        cds.append(c.drop(columns=["bakta_xrefs", "enhanced_xrefs"]))  # long GO lists, only needed for per_db
        clusters.append(cl.assign(genome=acc))
    if skipped:
        print(f"  skipped {len(skipped)} genomes with missing or ambiguous GFF3s: {', '.join(skipped)}")
    if not analysed:
        sys.exit(f"no genome has both GFF3s; is the drive mounted? ({config.ENHANCED_GFF_DIR})")
    cds, clusters = pd.concat(cds, ignore_index=True), pd.concat(clusters, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)

    write(cds, "bgc_cds.tsv")
    write(clusters[["genome", *clusters.columns.drop("genome")]], "clusters.tsv")
    pg = bgc.per_genome(cds, clusters, analysed)
    write(pg, "per_genome.tsv")
    write(pd.concat(dbs, ignore_index=True), "dbxref_by_db.tsv")
    write(bgc.summarise(pg), "summary.tsv")
    write(crosscheck(clusters, analysed), "crosscheck_cluster_counts.tsv")


def crosscheck(clusters: pd.DataFrame, genomes) -> pd.DataFrame:
    """antiSMASH clusters parsed here vs the pipeline summary."""
    n = clusters.groupby("genome").size().reindex(genomes, fill_value=0).rename("parsed_antiSMASH")
    summary = genome_list().set_index("File")[["Enhanced_BGCs"]]
    return n.to_frame().join(summary).rename_axis("genome").reset_index()


EXPORTS = {
    "summary.tsv": "summary.csv",
    "per_genome.tsv": "per_genome.csv",
    "dbxref_by_db.tsv": "dbxref_by_db.csv",
    "bgc_cds.tsv": "bgc_cds.csv",
    "clusters.tsv": "clusters.csv",
}


def export():
    config.EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    for src, dst in EXPORTS.items():
        df = pd.read_csv(OUT / src, sep="\t", keep_default_na=False)
        df.to_csv(config.EXPORT_DIR / dst, index=False)
        print(f"  exported {dst} ({len(df)} rows)")


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    steps = {"build": [build], "export": [export], "all": [build, export]}
    if step not in steps:
        sys.exit(f"usage: {sys.argv[0]} {{{','.join(steps)}}}")
    for f in steps[step]:
        f()
