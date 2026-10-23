"""KEGG analysis workflow.

Run from the repo root:
  python kegg_analyses/run_kegg_eval.py {validate,build,reconstruct,trace,review,export,table1,all}
  python kegg_analyses/run_kegg_eval.py table2     # after curating data/raw/kegg_analyses/review_candidates.csv
Requires the reference snapshot from fetch_kegg_reference.py. Working tables go to kegg_analyses/output/;
export copies the manuscript subset as CSV to data/raw/kegg_analyses/; Tables 1 and 2 go to tables/.
"""
import hashlib
import sys

import pandas as pd

import config
import metrics
import profiles
import review as rv
import traceability as tr
from modules import evaluate_panel, kos_in, parse

OUT = config.OUT_DIR


def write(df: pd.DataFrame, name: str):
    df.to_csv(OUT / name, sep="\t", index=False)
    print(f"  wrote {name} ({len(df)} rows)")


def validate() -> pd.DataFrame:
    m = profiles.load_manifest()
    errors, prov = [], []
    for g in m.itertuples():
        if not g.bakta_gff or not g.enhanced_gff:
            errors.append(f"{g.accession}: missing or ambiguous GFF3 (bakta={bool(g.bakta_gff)}, enhanced={bool(g.enhanced_gff)})")
            continue
        bakta, enh = profiles.read_cds(g.bakta_gff), profiles.read_cds(g.enhanced_gff)
        orphans = bakta.keys() - enh.keys()
        dropped = sum(bool(set(profiles.KO.findall(b.get("Dbxref", ""))) - set(profiles.KO.findall(enh[i].get("Dbxref", ""))))
                      for i, b in bakta.items() if i in enh)
        if orphans:
            errors.append(f"{g.accession}: {len(orphans)} Bakta CDS missing from enhanced GFF3, e.g. {sorted(orphans)[:3]}")
        if dropped:
            errors.append(f"{g.accession}: {dropped} proteins lose a Bakta KO in the enhanced GFF3 (A ⊆ D assumption broken)")
        for tool, path, n in [("bakta", g.bakta_gff, len(bakta)), ("enhanced", g.enhanced_gff, len(enh))]:
            header = " | ".join(x.strip() for x in open(path) if x.startswith("#") and not x.startswith("##sequence"))
            prov.append({"genome": g.accession, "file": tool, "path": path, "n_cds": n,
                         "md5": hashlib.md5(open(path, "rb").read()).hexdigest(), "header": header})
        prov.append({"genome": g.accession, "file": "enhanced_only_cds", "n_cds": len(enh.keys() - bakta.keys())})
    OUT.mkdir(parents=True, exist_ok=True)
    write(m, "genomes.tsv")
    write(pd.DataFrame(prov), "provenance.tsv")
    print(f"validate: {len(m)} genomes, {len(errors)} errors")
    if errors:
        sys.exit("\n".join(errors))
    return m


def load_trees() -> dict:
    d = pd.read_csv(config.KEGG_DIR / "module_definitions.tsv", sep="\t")
    return {r.module: parse(r.definition) for r in d.itertuples()}


def tiers(trees: dict) -> dict:
    return {"T1": set().union(*map(kos_in, trees.values())),
            "T2": set(pd.read_csv(config.KEGG_DIR / "metabolic_kos.tsv", sep="\t").ko)}


def build():
    m = pd.read_csv(OUT / "genomes.tsv", sep="\t")
    pfam2ipr = profiles.load_pfam2ipr()
    proteins = pd.concat([profiles.protein_table(g.accession, g.bakta_gff, g.enhanced_gff, pfam2ipr)
                          for g in m.itertuples()], ignore_index=True)
    pko = profiles.protein_ko(proteins)
    write(proteins, "proteins.tsv")
    write(pko, "protein_ko.tsv")

    t = tiers(load_trees())
    long, wide = metrics.coverage(proteins, pko, t)
    recovery = metrics.hypothetical_recovery(proteins, t)
    write(long, "coverage_per_state.tsv")
    write(wide, "coverage_per_genome.tsv")
    write(recovery, "hypothetical_recovery.tsv")
    write(metrics.summarise(wide, recovery), "table1_global.tsv")

    # Cross-check: same order of magnitude as the pipeline's Dbxref counts (occurrences, not distinct KOs)
    summary = pd.read_csv(config.PIPELINE_SUMMARY).set_index("File")
    distinct = pd.DataFrame({"distinct_A": pko[pko.in_bakta].groupby("genome").ko.nunique(),
                             "distinct_D": pko.groupby("genome").ko.nunique()})
    check = distinct.join(summary[["Initial_KEGG_entries", "Enhanced_KEGG_entries"]], how="left")
    check["ratio_A"] = check.Initial_KEGG_entries / check.distinct_A
    check["ratio_D"] = check.Enhanced_KEGG_entries / check.distinct_D
    write(check.reset_index(names="genome"), "crosscheck_pipeline_summary.tsv")
    print("IPS adds no KOs: state B = A and D = C for all KO metrics.")


def reconstruct():
    m = pd.read_csv(OUT / "genomes.tsv", sep="\t")
    pko = pd.read_csv(OUT / "protein_ko.tsv", sep="\t")
    trees = load_trees()
    write(metrics.reconstruct(pko, m.accession.tolist(), trees), "modules_per_genome.tsv")

    # Sanity check: the evaluator on KEGG's own organism annotations vs KEGG's module assignments
    rows = []
    for org in config.BENCHMARK_ORGS.values():
        kos = set((config.KEGG_DIR / f"{org}_kos.txt").read_text().split())
        listed = set((config.KEGG_DIR / f"{org}_modules.txt").read_text().split())
        for module, r in evaluate_panel(trees, kos).items():
            rows.append({"org": org, "module": module, "status": r["status"], "kegg_lists_module": module in listed,
                         "agree": r["complete"] == (module in listed) or module in config.COMPOSITE})
    check = pd.DataFrame(rows)
    write(check, "kegg_reference_check.tsv")
    if not check.agree.all():
        print("WARNING: evaluator disagrees with KEGG for\n" + check[~check.agree].to_string(index=False))


def module_names() -> dict:
    d = pd.read_csv(config.KEGG_DIR / "module_definitions.tsv", sep="\t")
    return dict(zip(d.module, d.name))


def ko_names() -> dict:
    d = pd.read_csv(config.KEGG_DIR / "ko_names.tsv", sep="\t")
    return dict(zip(d.ko, d.name))


def load(name: str) -> pd.DataFrame:
    return pd.read_csv(OUT / name, sep="\t", keep_default_na=False)


def trace():
    genomes, pko, proteins = load("genomes.tsv"), load("protein_ko.tsv"), load("proteins.tsv")
    trees = load_trees()
    supp, contributions = tr.trace(genomes, pko, proteins, trees, module_names())
    write(supp, "supp_genome_module.tsv")
    write(contributions, "module_contributions.tsv")
    write(tr.nitrogen_table(genomes, pko, proteins, supp), "nitrogen_table.tsv")
    write(tr.benchmark(pko, supp, trees), "benchmark_syz_cyu.tsv")
    write(tr.tca_new(pko, proteins, ko_names()), "tca_new_assignments.tsv")
    print(f"  {int(supp.requires_review.sum())} genome–module rows flagged for review")


CANDIDATES = config.EXPORT_DIR / "review_candidates.csv"


def review():
    scored = rv.score_pairs(load("protein_ko.tsv"), load("proteins.tsv"), load("module_contributions.tsv"),
                            load("genomes.tsv"), ko_names())
    table = rv.candidates_table(rv.select(scored))
    if CANDIDATES.exists():  # never overwrite manual curation
        table.to_csv(OUT / "review_candidates_new.csv", index=False)
        print(f"  kept curated {CANDIDATES.name}; fresh list written to output/review_candidates_new.csv")
    else:
        config.EXPORT_DIR.mkdir(parents=True, exist_ok=True)
        table.to_csv(CANDIDATES, index=False)
        print(f"  wrote {CANDIDATES.name} ({len(table)} candidates) for curation")


def table1():
    t = metrics.table1_manuscript(pd.read_csv(OUT / "table1_global.tsv", sep="\t"))
    config.TABLES_DIR.mkdir(exist_ok=True)
    t.to_csv(config.TABLES_DIR / "kegg_table1.csv", index=False)
    print(f"  wrote tables/kegg_table1.csv ({len(t)} rows)")


def table2():
    t = rv.table2(pd.read_csv(CANDIDATES, keep_default_na=False, dtype=str))
    if not 5 <= len(t) <= 10:
        print(f"  note: {len(t)} rows selected; expected 5–10")
    config.TABLES_DIR.mkdir(exist_ok=True)
    t.to_csv(config.TABLES_DIR / "kegg_table2.csv", index=False)
    print(f"  wrote tables/kegg_table2.csv ({len(t)} rows)")


# output table -> manuscript CSV (supplementary material, data availability, plot inputs)
EXPORTS = {
    "genomes.tsv": "genomes.csv",
    "table1_global.tsv": "table1_global.csv",
    "coverage_per_genome.tsv": "coverage_per_genome.csv",
    "hypothetical_recovery.tsv": "hypothetical_recovery.csv",
    "modules_per_genome.tsv": "modules_per_genome.csv",
    "protein_ko.tsv": "protein_ko.csv",
    "supp_genome_module.tsv": "supp_genome_module.csv",
    "nitrogen_table.tsv": "nitrogen_table.csv",
    "benchmark_syz_cyu.tsv": "benchmark_syz_cyu.csv",
}


def export():
    config.EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    for src, dst in EXPORTS.items():
        df = pd.read_csv(OUT / src, sep="\t", keep_default_na=False)
        if src == "genomes.tsv":
            df = df.drop(columns=["bakta_gff", "enhanced_gff"])  # local paths, meaningless outside this machine
        df.to_csv(config.EXPORT_DIR / dst, index=False)
        print(f"  exported {dst} ({len(df)} rows)")
    t2 = tier2_maps()
    t2.to_csv(config.EXPORT_DIR / "tier2_maps.csv", index=False)
    print(f"  exported tier2_maps.csv ({len(t2)} rows)")


def tier2_maps() -> pd.DataFrame:
    """T2 definition per map, with the distinct KOs observed across all genomes in states A and D."""
    m = pd.read_csv(config.KEGG_DIR / "metabolic_kos.tsv", sep="\t")
    names = pd.read_csv(config.KEGG_DIR / "map_names.tsv", sep="\t")
    pko = pd.read_csv(OUT / "protein_ko.tsv", sep="\t")
    obs = m.assign(in_bakta=m.ko.isin(set(pko.ko[pko.in_bakta])), in_integrated=m.ko.isin(set(pko.ko)))
    t = obs.groupby("map").agg(kos_in_map=("ko", "size"), kos_bakta=("in_bakta", "sum"),
                               kos_integrated=("in_integrated", "sum")).reset_index()
    t["kos_added"] = t.kos_integrated - t.kos_bakta
    return names.merge(t, on="map", how="right")[["map", "name", *t.columns[1:]]]


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    steps = {"validate": [validate], "build": [build], "reconstruct": [reconstruct], "trace": [trace],
             "review": [review], "export": [export], "table1": [table1], "table2": [table2],
             "all": [validate, build, reconstruct, trace, review, export, table1]}
    if step not in steps:
        sys.exit(f"usage: {sys.argv[0]} {{{','.join(steps)}}}")
    for f in steps[step]:
        f()
