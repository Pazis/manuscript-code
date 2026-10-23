"""Snapshot the KEGG and InterPro data the analysis needs into kegg_analyses/reference/.

Downstream steps read only the snapshot. Rerunning this overwrites it with the current releases and can change the results.
"""
import datetime
import gzip
import time
import urllib.request
import xml.etree.ElementTree as ET

from config import KEGG_DIR, PANEL, BENCHMARK_ORGS, PFAM2IPR, is_metabolic_map

REST = "https://rest.kegg.jp"
INTERPRO_XML = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz"


def get(path: str) -> str:
    time.sleep(0.4)  # KEGG asks for <3 requests/s
    with urllib.request.urlopen(f"{REST}/{path}", timeout=120) as r:
        return r.read().decode()


def definition(entry: str) -> tuple[str, str]:
    name, lines, inside = "", [], False
    for line in entry.splitlines():
        key = line[:12].strip()
        if key == "NAME":
            name = line[12:].strip()
        if key == "DEFINITION":
            inside = True
        elif key:
            inside = False
        if inside:
            lines.append(line[12:].strip())
    return name, " ".join(lines)


def fetch_pfam2interpro():
    """Pfam accession/short name -> integrating InterPro entry, streamed from interpro.xml.gz."""
    rows, release = [], ""
    with urllib.request.urlopen(INTERPRO_XML, timeout=600) as r, gzip.open(r) as xml:
        for _, el in ET.iterparse(xml):
            if el.tag == "dbinfo" and el.get("dbname") == "INTERPRO":
                release = f"InterPro {el.get('version')} ({el.get('file_date')})"
            elif el.tag == "interpro":
                rows += [(x.get("dbkey"), x.get("name"), el.get("id"))
                         for x in el.iterfind("member_list/db_xref") if x.get("db") == "PFAM"]
                el.clear()
    PFAM2IPR.parent.mkdir(parents=True, exist_ok=True)
    PFAM2IPR.write_text("pfam\tpfam_name\tinterpro\n" + "".join(f"{a}\t{n}\t{i}\n" for a, n, i in sorted(rows)))
    (PFAM2IPR.parent / "RELEASE.txt").write_text(f"fetched {datetime.date.today()}\n{release}\n")


def main():
    fetch_pfam2interpro()
    raw = KEGG_DIR / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    rows = []
    for m in PANEL:
        text = get(f"get/{m}")
        (raw / f"{m}.txt").write_text(text)
        rows.append((m, *definition(text)))
    (KEGG_DIR / "module_definitions.tsv").write_text(
        "module\tname\tdefinition\n" + "".join(f"{m}\t{n}\t{d}\n" for m, n, d in rows))

    (KEGG_DIR / "ko_names.tsv").write_text("ko\tname\n" + get("list/ko").replace("ko:", ""))
    (KEGG_DIR / "map_names.tsv").write_text("map\tname\n" + get("list/pathway"))

    pairs = set()
    for line in get("link/pathway/ko").splitlines():
        ko, path = line.split("\t")
        path = path.removeprefix("path:")
        if path.startswith("map") and is_metabolic_map(path):
            pairs.add((ko.removeprefix("ko:"), path))
    (KEGG_DIR / "metabolic_kos.tsv").write_text(
        "ko\tmap\n" + "".join(f"{k}\t{p}\n" for k, p in sorted(pairs)))

    for org in BENCHMARK_ORGS.values():
        kos = sorted({l.split("\t")[1].removeprefix("ko:") for l in get(f"link/ko/{org}").splitlines()})
        (KEGG_DIR / f"{org}_kos.txt").write_text("\n".join(kos) + "\n")
        mods = sorted({l.split("\t")[1].split("_")[-1] for l in get(f"link/module/{org}").splitlines()})
        (KEGG_DIR / f"{org}_modules.txt").write_text("\n".join(mods) + "\n")

    (KEGG_DIR / "RELEASE.txt").write_text(
        f"fetched {datetime.date.today()}\n\n" + get("info/module") + "\n" + get("info/ko"))


if __name__ == "__main__":
    main()
