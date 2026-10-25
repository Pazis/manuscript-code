"""GFF3 parsing and per-CDS annotation gain inside antiSMASH clusters.

Reannota appends antiSMASH (and GECCO) clusters to the enhanced GFF3 as `biosynthetic-gene-cluster`
features (source `antiSMASH`, `antiSMASH:<version>` or `GECCO`), each followed by copies of its
ORFs as extra CDS rows. Only antiSMASH clusters are analysed. Membership is taken from cluster
coordinates, not from those copies: a Bakta CDS is in a cluster when it lies entirely inside it on the same contig.
"""
from urllib.parse import unquote

import pandas as pd

import config

BGC_SOURCES = ("antiSMASH", "GECCO")  # their CDS rows are copies of real CDS and are skipped


def read_gff(path) -> tuple[dict, list]:
    """(CDS by ID, antiSMASH clusters)."""
    cds, clusters = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            col = line.rstrip("\n").split("\t")
            if len(col) != 9:
                continue
            attrs = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in col[8].split(";") if "=" in kv)}
            from_bgc_tool = col[1].startswith(BGC_SOURCES)
            if col[2] == "biosynthetic-gene-cluster" and col[1].startswith("antiSMASH"):
                clusters.append({"bgc": attrs.get("ID", ""), "seqid": col[0], "start": int(col[3]),
                                 "end": int(col[4]), "bgc_product": attrs.get("product", "")})
            elif col[2] == "CDS" and not from_bgc_tool:
                if attrs["ID"] in cds:
                    raise ValueError(f"duplicate CDS ID {attrs['ID']} in {path}")
                cds[attrs["ID"]] = {"seqid": col[0], "start": int(col[3]), "end": int(col[4]), **attrs}
    return cds, clusters


COUNTED_DBS = {"KEGG", "GO", "COG", "EC"}


def dbxrefs(attrs: dict) -> set:
    """KEGG, GO, COG and EC entries only. COG category letters (Bakta `COG:L`, eggNOG `COG_cat:L`) are dropped."""
    out = set()
    for x in filter(None, attrs.get("Dbxref", "").split(",")):
        x = x.replace("EC:EC:", "EC:", 1)  # enhanced GFF3s sometimes double the EC prefix
        db, _, acc = x.partition(":")
        if db in COUNTED_DBS and (db != "COG" or acc.startswith("COG")):
            out.add(x)
    return out


def db_of(xref: str) -> str:
    return xref.split(":", 1)[0]


COLUMNS = ["genome", "protein", "seqid", "start", "end", "bgcs", "bgc_products", "bakta_product", "enhanced_product",
           "is_hypothetical", "recovered", "dbxrefs_bakta", "dbxrefs_enhanced", "dbxrefs_added", "dbxrefs_lost",
           "added_dbs", "improved", "bakta_xrefs", "enhanced_xrefs"]


def cds_table(genome: str, bakta_gff, enhanced_gff) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per Bakta CDS inside an antiSMASH cluster, with its annotation gain; plus the genome's clusters."""
    bakta, _ = read_gff(bakta_gff)
    enh, clusters = read_gff(enhanced_gff)
    missing = bakta.keys() - enh.keys()
    if missing:
        raise ValueError(f"{genome}: {len(missing)} Bakta CDS missing from enhanced GFF3, e.g. {sorted(missing)[:3]}")
    unknown = {c["seqid"] for c in clusters} - {b["seqid"] for b in bakta.values()}
    if unknown:
        raise ValueError(f"{genome}: clusters on contigs with no Bakta CDS: {sorted(unknown)[:3]}")

    by_contig = {}
    for c in clusters:
        by_contig.setdefault(c["seqid"], []).append(c)
    rows = []
    for pid, b in bakta.items():
        inside = [c for c in by_contig.get(b["seqid"], []) if c["start"] <= b["start"] and b["end"] <= c["end"]]
        if not inside:
            continue
        e = enh[pid]
        bx, ex = dbxrefs(b), dbxrefs(e)
        added = ex - bx
        hyp = b.get("product", "") == config.HYPOTHETICAL
        recovered = hyp and e.get("product", "") not in config.UNINFORMATIVE
        rows.append({"genome": genome, "protein": pid, "seqid": b["seqid"], "start": b["start"], "end": b["end"],
                     "bgcs": ",".join(c["bgc"] for c in inside),  # overlapping protoclusters: CDS still counted once
                     "bgc_products": ",".join(c["bgc_product"] for c in inside),
                     "bakta_product": b.get("product", ""), "enhanced_product": e.get("product", ""),
                     "is_hypothetical": hyp, "recovered": recovered,
                     "dbxrefs_bakta": len(bx), "dbxrefs_enhanced": len(ex),
                     "dbxrefs_added": len(added), "dbxrefs_lost": len(bx - ex),
                     "added_dbs": ",".join(sorted({db_of(x) for x in added})),
                     "improved": recovered or bool(added),
                     "bakta_xrefs": ",".join(sorted(bx)), "enhanced_xrefs": ",".join(sorted(ex))})  # for per_db
    return (pd.DataFrame(rows, columns=COLUMNS),
            pd.DataFrame(clusters, columns=["bgc", "seqid", "start", "end", "bgc_product"]))


def per_genome(cds: pd.DataFrame, clusters: pd.DataFrame, genomes) -> pd.DataFrame:
    """Counts per genome; genomes without antiSMASH clusters get zeros. Hypothetical counts are in CDS space."""
    g = cds.groupby("genome").agg(
        n_cds=("protein", "size"), n_hypothetical=("is_hypothetical", "sum"), n_recovered=("recovered", "sum"),
        n_gained_dbxref=("dbxrefs_added", lambda s: int((s > 0).sum())), n_improved=("improved", "sum"),
        dbxrefs_bakta=("dbxrefs_bakta", "sum"), dbxrefs_enhanced=("dbxrefs_enhanced", "sum"),
        dbxrefs_added=("dbxrefs_added", "sum"), dbxrefs_lost=("dbxrefs_lost", "sum"))
    g.insert(0, "n_clusters", clusters.groupby("genome").size())
    t = g.reindex(genomes, fill_value=0).fillna(0).astype(int).rename_axis("genome").reset_index()
    t["pct_hyp_recovered"] = 100 * t.n_recovered / t.n_hypothetical.where(t.n_hypothetical > 0)
    t["pct_cds_improved"] = 100 * t.n_improved / t.n_cds.where(t.n_cds > 0)
    return t


def per_db(cds: pd.DataFrame) -> pd.DataFrame:
    """Dbxref entries per genome and database (GO, InterPro, KEGG, COG, EC, ...)."""
    out = []
    for state, col in [("bakta", "bakta_xrefs"), ("enhanced", "enhanced_xrefs")]:
        x = cds[["genome"]].assign(db=cds[col].str.split(",")).explode("db")
        x = x[x.db.notna() & (x.db != "")]
        out.append(x.assign(db=x.db.map(db_of), state=state).groupby(["genome", "db", "state"]).size().rename("n"))
    t = pd.concat(out).unstack("state", fill_value=0).reindex(columns=["bakta", "enhanced"], fill_value=0).reset_index()
    t.columns.name = None
    t["added"] = t.enhanced - t.bakta
    return t


def summarise(pg: pd.DataFrame) -> pd.DataFrame:
    """Totals across genomes, pooled %, and median [IQR] of the per-genome percentages (one row)."""
    r = {"n_genomes": len(pg), "n_genomes_with_clusters": int((pg.n_clusters > 0).sum())}
    for c in ["n_clusters", "n_cds", "n_hypothetical", "n_recovered", "n_gained_dbxref", "n_improved",
              "dbxrefs_bakta", "dbxrefs_enhanced", "dbxrefs_added", "dbxrefs_lost"]:
        r[c] = int(pg[c].sum())
    r["pct_hyp_recovered_pooled"] = 100 * r["n_recovered"] / r["n_hypothetical"] if r["n_hypothetical"] else float("nan")
    r["pct_cds_improved_pooled"] = 100 * r["n_improved"] / r["n_cds"] if r["n_cds"] else float("nan")
    for c in ["pct_hyp_recovered", "pct_cds_improved"]:
        r[f"{c}_median"], r[f"{c}_q1"], r[f"{c}_q3"] = pg[c].median(), pg[c].quantile(.25), pg[c].quantile(.75)
    return pd.DataFrame([r])
