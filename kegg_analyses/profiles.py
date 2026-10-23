"""Manifest, GFF3 parsing and protein–KO evidence with provenance.

Provenance is derived from two GFF3 files per genome that share CDS IDs:
  Bakta GFF3     -> baseline (state A). KOs are KEGG: Dbxrefs.
  enhanced GFF3  -> Reannota merge of Bakta + InterProScan + eggNOG-mapper (state D).
InterProScan assigns no KOs, so a KO on the enhanced CDS that Bakta did not give it came from eggNOG-mapper.
"""
import re
from urllib.parse import unquote

import pandas as pd

import config

KO = re.compile(r"KEGG:(K\d{5})")
EC = re.compile(r"EC:(?:EC:)?(\d+\.[\d\-]+\.[\d\-]+\.n?[\d\-]+)")
IPR = re.compile(r"InterPro:(IPR\d{6})")
EGGNOG_PFAM = "Eggnog PFAM comment:"
# Products that replace "hypothetical protein" without naming a function (MobiDB-lite output)
UNINFORMATIVE = {config.HYPOTHETICAL, "consensus disorder prediction"}


def load_manifest() -> pd.DataFrame:
    m = pd.read_excel(config.GENOMES_XLSX, sheet_name=config.MANIFEST_SHEET)
    m.columns = ["group", "accession", "organism", "family", "assembly_level", "role"]
    m["substitution_note"] = ""
    genomes = pd.read_excel(config.GENOMES_XLSX, sheet_name="Cyanobacterial Genomes").set_index("accession")
    for old, (new, organism, note) in config.SUBSTITUTIONS.items():
        i = m.index[m.accession == old]
        if len(i) != 1:
            raise ValueError(f"substitution target {old} not found once in {config.MANIFEST_SHEET}")
        m.loc[i, ["accession", "organism", "substitution_note"]] = [new, organism, note]
        m.loc[i, ["family", "assembly_level"]] = genomes.loc[new, ["family", "assembly_level"]].values
    if m.accession.duplicated().any():
        raise ValueError(f"duplicate accessions: {m.accession[m.accession.duplicated()].tolist()}")
    m["kegg_org"] = m.accession.map(config.BENCHMARK_ORGS).fillna("")
    m["bakta_gff"] = m.accession.map(lambda a: find_one(config.BAKTA_GFF_DIR, a))
    m["enhanced_gff"] = m.accession.map(lambda a: find_one(config.ENHANCED_GFF_DIR, a))
    return m


def find_one(directory, accession) -> str:
    hits = sorted(directory.glob(f"{accession}_*.gff3"))
    return str(hits[0]) if len(hits) == 1 else ""  # missing or ambiguous: reported by validate


def read_cds(path) -> dict:
    """CDS attributes by ID. antiSMASH re-emits BGC genes as extra CDS features (ID=as:...); those are skipped."""
    cds = {}
    with open(path) as fh:
        for line in fh:
            col = line.rstrip("\n").split("\t")
            if len(col) != 9 or col[2] != "CDS" or col[1] == "antiSMASH":
                continue
            attrs = dict(kv.split("=", 1) for kv in col[8].split(";") if "=" in kv)
            if attrs["ID"] in cds:
                raise ValueError(f"duplicate CDS ID {attrs['ID']} in {path}")
            cds[attrs["ID"]] = {k: unquote(v) for k, v in attrs.items()}
    return cds


def eggnog_pfams(note: str) -> set:
    # Reannota writes "Eggnog PFAM comment:<pfam>,<pfam>,..." as the last Note item
    _, _, tail = note.partition(EGGNOG_PFAM)
    return {p for p in tail.split(",") if p}


def protein_table(genome: str, bakta_gff, enhanced_gff, pfam2ipr: dict) -> pd.DataFrame:
    bakta, enh = read_cds(bakta_gff), read_cds(enhanced_gff)
    rows = []
    for pid, b in bakta.items():
        e = enh[pid]
        b_x, e_x, note = b.get("Dbxref", ""), e.get("Dbxref", ""), e.get("Note", "")
        bakta_kos, enh_kos = set(KO.findall(b_x)), set(KO.findall(e_x))
        iprs = set(IPR.findall(e_x))
        pfams = eggnog_pfams(note)
        rows.append({
            "genome": genome, "protein": pid,
            "bakta_product": b.get("product", ""), "bakta_gene": b.get("gene", ""),
            "enhanced_product": e.get("product", ""), "enhanced_name": e.get("Name", ""),
            "is_hypothetical": b.get("product", "") == config.HYPOTHETICAL,
            "bakta_kos": ",".join(sorted(bakta_kos)), "enhanced_kos": ",".join(sorted(enh_kos)),
            "new_kos": ",".join(sorted(enh_kos - bakta_kos)), "dropped_kos": ",".join(sorted(bakta_kos - enh_kos)),
            "bakta_ec": ",".join(sorted(set(EC.findall(b_x)))), "enhanced_ec": ",".join(sorted(set(EC.findall(e_x)))),
            "interpro": ",".join(sorted(iprs)), "eggnog_pfams": ",".join(sorted(pfams)),
            "has_ips": bool(iprs) or "InterPro detected" in note,
            "has_eggnog": bool(pfams) or bool(enh_kos - bakta_kos),
            # eggNOG's Pfams overlap the InterProScan hits (compared via InterPro integration)
            "ips_confirmed": bool({pfam2ipr.get(p) for p in pfams} & iprs),
            "named_after": e.get("product", "") not in UNINFORMATIVE,
        })
    return pd.DataFrame(rows)


def protein_ko(proteins: pd.DataFrame) -> pd.DataFrame:
    """One row per (genome, protein, KO) in the integrated profile, with provenance."""
    rows = []
    for p in proteins.itertuples():
        bakta = set(filter(None, p.bakta_kos.split(",")))
        for ko in filter(None, p.enhanced_kos.split(",")):
            rows.append({"genome": p.genome, "protein": p.protein, "ko": ko,
                         "in_bakta": ko in bakta, "in_eggnog": ko not in bakta,
                         "provenance": "Bakta" if ko in bakta else "eggNOG",
                         "ips_confirmed": p.ips_confirmed, "is_hypothetical": p.is_hypothetical,
                         "conflict": ko not in bakta and bool(bakta)})
    return pd.DataFrame(rows)


def load_pfam2ipr() -> dict:
    t = pd.read_csv(config.PFAM2IPR, sep="\t")
    return dict(zip(t.pfam_name, t.interpro))
