"""antiSMASH vs GECCO: cluster consolidation and boundary overlap per genome (Table S8, Figures 6-7).

Reads the Reannota enhanced GFF3s (`<genome>/Reannota_New/enhanced/<genome>_enhanced.gff3`), which
end with one `biosynthetic-gene-cluster` row per antiSMASH region (source `antiSMASH`) and per GECCO
cluster (source `GECCO`). Everything is computed from those coordinates.

Run from the repo root:
  python bgc_analyses/bgc_consolidation.py            # tables + update data/processed/BGC_analyses_tables.xlsx
  python bgc_analyses/bgc_consolidation.py --no-workbook
  python -m pytest bgc_analyses/tests

Unit of comparison: one antiSMASH **region** vs one GECCO cluster. GFF3s written before the ReAnnota
fix (antiSMASH protoclusters, IDs `as:<contig>_bgc<n>`) are refused, because protoclusters overlap each
other and inflate every count and overlap below.

Per genome (columns of sheet "S8 BGC consolidation"; the first eleven keep their old names so the
figure scripts still work):
  AntiSMASH Count / GECCO Count   clusters reported by each tool
  Overlapping Clusters            antiSMASH regions that overlap >= 1 GECCO cluster (>= 1 bp, same contig)
  GECCO Products                  GECCO cluster types, most frequent first
  Total Overlap (bp)              length of (union of antiSMASH regions) intersect (union of GECCO clusters)
  % AntiSMASH Covered             Total Overlap / union length of antiSMASH regions x 100
  % GECCO Covered                 Total Overlap / union length of GECCO clusters x 100
  Inter-tool Consensus            Overlapping Clusters / AntiSMASH Count
  Avg Overlap Efficiency          mean Jaccard (shared bp / union bp) over overlapping region-cluster pairs; 0..1
  Consolidation Ratio             AntiSMASH Count / GECCO Count
  AntiSMASH Locus (bp)            union length of antiSMASH regions (denominator of % AntiSMASH Covered)
  GECCO Locus (bp)                union length of GECCO clusters (denominator of % GECCO Covered)
  GECCO Matched                   GECCO clusters that overlap >= 1 antiSMASH region
  GECCO Consensus                 GECCO Matched / GECCO Count
  AntiSMASH Only / GECCO Only     clusters with no counterpart in the other tool
  Regions with >1 GECCO           antiSMASH regions overlapping two or more GECCO clusters (GECCO splits)
  GECCO with >1 Region            GECCO clusters overlapping two or more antiSMASH regions (antiSMASH splits)
  Matched Consolidation Ratio     Overlapping Clusters / GECCO Matched: the ratio among co-detected clusters,
                                  i.e. without clusters that only one tool found
"""
import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATASET_DIR = Path("/media/1GBDrive/cyano-annotator/Dataset_Good_Quality")
GFF_GLOB = "*/Reannota_New/enhanced/*_enhanced.gff3"
WORKBOOK = ROOT / "data/processed/BGC_analyses_tables.xlsx"
SHEET = "S8 BGC consolidation"
SUMMARY_SHEET = "S8 summary"
OUT_DIR = HERE / "output"

OLD_COLUMNS = ["Locus Tag", "AntiSMASH Count", "GECCO Count", "Overlapping Clusters", "GECCO Products",
               "Total Overlap (bp)", "% AntiSMASH Covered", "% GECCO Covered", "Inter-tool Consensus",
               "Avg Overlap Efficiency", "Consolidation Ratio"]
NEW_COLUMNS = ["AntiSMASH Locus (bp)", "GECCO Locus (bp)", "GECCO Matched", "GECCO Consensus", "AntiSMASH Only", "GECCO Only", "Regions with >1 GECCO",
               "GECCO with >1 Region", "Matched Consolidation Ratio"]
COLUMNS = OLD_COLUMNS + NEW_COLUMNS

PROTOCLUSTER_ID = re.compile(r"^as:.+_bgc\d+$")


# --------------------------------------------------------------------------- parsing
def _attr(col9: str, key: str) -> str:
    """Attribute value. Reannota writes multi-type products unescaped (product=NRP;Polyketide;tool=GECCO),
    so `product` is read up to `;tool=`."""
    if key == "product":
        m = re.search(r"(?:^|;)product=(.*?)(?:;tool=|$)", col9)
    else:
        m = re.search(rf"(?:^|;){key}=([^;]*)", col9)
    return m.group(1) if m else ""


def read_clusters(gff_path) -> tuple[list, list, set]:
    """(antiSMASH regions, GECCO clusters, contigs carrying genome features) from one enhanced GFF3."""
    regions, gecco, contigs = [], [], set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) != 9:
                continue
            if col[2] != "biosynthetic-gene-cluster":
                if not col[1].startswith(("antiSMASH", "GECCO")):
                    contigs.add(col[0])
                continue
            rec = {"id": _attr(col[8], "ID"), "seqid": col[0], "start": int(col[3]), "end": int(col[4]),
                   "product": _attr(col[8], "product")}
            if rec["start"] > rec["end"]:
                raise ValueError(f"{gff_path}: cluster {rec['id']} has start > end")
            if col[1].startswith("antiSMASH"):
                if PROTOCLUSTER_ID.match(rec["id"]):
                    raise ValueError(f"{gff_path}: antiSMASH protocluster IDs ({rec['id']}); this GFF3 predates the "
                                     "ReAnnota region fix. Rerun ReAnnota and use the Reannota_New outputs.")
                regions.append(rec)
            elif col[1].startswith("GECCO"):
                gecco.append(rec)
    return regions, gecco, contigs


# --------------------------------------------------------------------------- interval arithmetic (1-based, inclusive)
def _merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def union_length(clusters) -> int:
    by_contig = defaultdict(list)
    for c in clusters:
        by_contig[c["seqid"]].append((c["start"], c["end"]))
    return sum(e - s + 1 for iv in by_contig.values() for s, e in _merge(iv))


def intersection_length(a, b) -> int:
    """Length of (union of a) intersect (union of b), contig by contig."""
    total = 0
    contigs = {c["seqid"] for c in a} & {c["seqid"] for c in b}
    for ctg in contigs:
        ua = _merge([(c["start"], c["end"]) for c in a if c["seqid"] == ctg])
        ub = _merge([(c["start"], c["end"]) for c in b if c["seqid"] == ctg])
        i = j = 0
        while i < len(ua) and j < len(ub):
            lo, hi = max(ua[i][0], ub[j][0]), min(ua[i][1], ub[j][1])
            if lo <= hi:
                total += hi - lo + 1
            if ua[i][1] < ub[j][1]:
                i += 1
            else:
                j += 1
    return total


def pair_overlap(a, b) -> int:
    if a["seqid"] != b["seqid"]:
        return 0
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)


# --------------------------------------------------------------------------- metrics
def genome_metrics(genome: str, regions: list, gecco: list) -> tuple[dict, list]:
    """One S8 row and the list of overlapping region-cluster pairs."""
    pairs = []
    for r in regions:
        for g in gecco:
            ov = pair_overlap(r, g)
            if ov:
                union = (r["end"] - r["start"] + 1) + (g["end"] - g["start"] + 1) - ov
                pairs.append({"genome": genome, "antismash_region": r["id"], "gecco_cluster": g["id"],
                              "seqid": r["seqid"], "region_start": r["start"], "region_end": r["end"],
                              "gecco_start": g["start"], "gecco_end": g["end"], "overlap_bp": ov,
                              "union_bp": union, "jaccard": ov / union,
                              "frac_region_covered": ov / (r["end"] - r["start"] + 1),
                              "frac_gecco_covered": ov / (g["end"] - g["start"] + 1),
                              "region_product": r["product"], "gecco_product": g["product"]})
    per_region = Counter(p["antismash_region"] for p in pairs)
    per_gecco = Counter(p["gecco_cluster"] for p in pairs)
    n_as, n_ge = len(regions), len(gecco)
    overlap = intersection_length(regions, gecco)
    u_as, u_ge = union_length(regions), union_length(gecco)
    types = Counter(t.strip() for g in gecco for t in g["product"].split(";") if t.strip())

    def ratio(a, b):
        return a / b if b else float("nan")

    row = {
        "Locus Tag": genome,
        "AntiSMASH Count": n_as,
        "GECCO Count": n_ge,
        "Overlapping Clusters": len(per_region),
        "GECCO Products": ", ".join(t for t, _ in types.most_common()),
        "Total Overlap (bp)": overlap,
        "% AntiSMASH Covered": 100 * ratio(overlap, u_as),
        "% GECCO Covered": 100 * ratio(overlap, u_ge),
        "Inter-tool Consensus": ratio(len(per_region), n_as),
        "Avg Overlap Efficiency": mean(p["jaccard"] for p in pairs) if pairs else float("nan"),
        "Consolidation Ratio": ratio(n_as, n_ge),
        "AntiSMASH Locus (bp)": u_as,
        "GECCO Locus (bp)": u_ge,
        "GECCO Matched": len(per_gecco),
        "GECCO Consensus": ratio(len(per_gecco), n_ge),
        "AntiSMASH Only": n_as - len(per_region),
        "GECCO Only": n_ge - len(per_gecco),
        "Regions with >1 GECCO": sum(1 for v in per_region.values() if v > 1),
        "GECCO with >1 Region": sum(1 for v in per_gecco.values() if v > 1),
        "Matched Consolidation Ratio": ratio(len(per_region), len(per_gecco)),
    }
    return row, pairs


def check_row(row: dict) -> list:
    """Sanity checks; returns problems (empty = fine)."""
    p = []
    for c in ["Inter-tool Consensus", "Avg Overlap Efficiency", "GECCO Consensus"]:
        v = row[c]
        if v == v and not 0 <= v <= 1 + 1e-9:
            p.append(f"{c}={v:.3f} outside 0..1")
    for c in ["% AntiSMASH Covered", "% GECCO Covered"]:
        v = row[c]
        if v == v and not 0 <= v <= 100 + 1e-9:
            p.append(f"{c}={v:.1f} outside 0..100")
    if row["Overlapping Clusters"] > row["AntiSMASH Count"] or row["GECCO Matched"] > row["GECCO Count"]:
        p.append("more matched clusters than clusters")
    if (row["Overlapping Clusters"] == 0) != (row["GECCO Matched"] == 0):
        p.append("overlap found from one side only")
    return p


def overlapping_within(clusters) -> int:
    """Pairs of clusters from the same tool that overlap each other (expected 0 for antiSMASH regions)."""
    n = 0
    for i, a in enumerate(clusters):
        for b in clusters[i + 1:]:
            n += bool(pair_overlap(a, b))
    return n


def summarise(df: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    """Cohort numbers for the manuscript. 'Plotted' = genomes with both tools and >= 1 overlap (Figures 6-7)."""
    plot = df[(df["AntiSMASH Count"] > 0) & (df["GECCO Count"] > 0) & (df["Overlapping Clusters"] > 0)]
    plot_pairs = pairs[pairs.genome.isin(plot["Locus Tag"])]
    s = lambda c: int(df[c].sum())  # noqa: E731
    rows = [
        ("Genomes analysed", len(df)),
        ("Genomes with antiSMASH regions", int((df["AntiSMASH Count"] > 0).sum())),
        ("Genomes with GECCO clusters", int((df["GECCO Count"] > 0).sum())),
        ("Genomes plotted (both tools, >= 1 overlap)", len(plot)),
        ("antiSMASH regions (all genomes)", s("AntiSMASH Count")),
        ("GECCO clusters (all genomes)", s("GECCO Count")),
        ("antiSMASH regions overlapping GECCO", s("Overlapping Clusters")),
        ("GECCO clusters overlapping antiSMASH", s("GECCO Matched")),
        ("antiSMASH-only regions", s("AntiSMASH Only")),
        ("GECCO-only clusters", s("GECCO Only")),
        ("antiSMASH regions with >1 GECCO cluster", s("Regions with >1 GECCO")),
        ("GECCO clusters spanning >1 antiSMASH region", s("GECCO with >1 Region")),
        ("Overlapping region-cluster pairs", len(pairs)),
        ("Inter-tool consensus, pooled (all genomes)", s("Overlapping Clusters") / s("AntiSMASH Count")),
        ("Inter-tool consensus, mean (plotted)", plot["Inter-tool Consensus"].mean()),
        ("Inter-tool consensus, median (plotted)", plot["Inter-tool Consensus"].median()),
        ("GECCO consensus, pooled (all genomes)", s("GECCO Matched") / s("GECCO Count")),
        ("Overlap efficiency, mean of per-genome means (plotted)", plot["Avg Overlap Efficiency"].mean()),
        ("Overlap efficiency, median over pairs (plotted)", plot_pairs.jaccard.median()),
        ("% antiSMASH locus covered by GECCO, median (plotted)", plot["% AntiSMASH Covered"].median()),
        ("% GECCO locus covered by antiSMASH, median (plotted)", plot["% GECCO Covered"].median()),
        ("% antiSMASH locus covered by GECCO, pooled bp (all genomes)",
         100 * s("Total Overlap (bp)") / s("AntiSMASH Locus (bp)")),
        ("% GECCO locus covered by antiSMASH, pooled bp (all genomes)",
         100 * s("Total Overlap (bp)") / s("GECCO Locus (bp)")),
        ("Consolidation ratio, mean (plotted)", plot["Consolidation Ratio"].mean()),
        ("Consolidation ratio, median (plotted)", plot["Consolidation Ratio"].median()),
        ("Consolidation ratio, pooled (all genomes)", s("AntiSMASH Count") / s("GECCO Count")),
        ("Matched consolidation ratio, mean (plotted)", plot["Matched Consolidation Ratio"].mean()),
        ("Matched consolidation ratio, pooled (all genomes)", s("Overlapping Clusters") / s("GECCO Matched")),
    ]
    return pd.DataFrame(rows, columns=["Metric", "Value"])


# --------------------------------------------------------------------------- driver
def accession(gff: Path) -> str:
    return "_".join(gff.name.split("_")[:2])


def build(dataset: Path):
    gffs = sorted(dataset.glob(GFF_GLOB))
    if not gffs:
        sys.exit(f"no enhanced GFF3s under {dataset}/{GFF_GLOB}; is the drive mounted?")
    rows, pairs, problems = [], [], []
    seen = Counter(accession(g) for g in gffs)
    for acc, n in seen.items():
        if n > 1:
            problems.append(f"{acc}: {n} GFF3s; all kept, check for duplicates")
    for gff in gffs:
        acc = accession(gff)
        regions, gecco, contigs = read_clusters(gff)
        for tool, cl in [("antiSMASH", regions), ("GECCO", gecco)]:
            stray = {c["seqid"] for c in cl} - contigs
            if stray:
                problems.append(f"{acc}: {tool} clusters on contigs absent from the annotation {sorted(stray)[:3]} "
                                "(contig names may differ between tools; overlaps would be missed)")
        if overlapping_within(regions):
            problems.append(f"{acc}: {overlapping_within(regions)} overlapping antiSMASH region pairs")
        row, pr = genome_metrics(acc, regions, gecco)
        problems += [f"{acc}: {p}" for p in check_row(row)]
        rows.append(row)
        pairs += pr
    df = pd.DataFrame(rows, columns=COLUMNS).sort_values("Locus Tag").reset_index(drop=True)
    pairs = pd.DataFrame(pairs)
    return df, pairs, summarise(df, pairs), problems


def rounded(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in ["% AntiSMASH Covered", "% GECCO Covered"]:
        out[c] = out[c].round(2)
    for c in ["Inter-tool Consensus", "Avg Overlap Efficiency", "Consolidation Ratio", "GECCO Consensus",
              "Matched Consolidation Ratio"]:
        out[c] = out[c].round(2)
    return out


def update_workbook(path: Path, df: pd.DataFrame, summary: pd.DataFrame):
    """Rewrite the S8 rows (keeping the header style), append the new columns, and (re)write the summary sheet.
    All other sheets are left untouched."""
    from copy import copy

    import openpyxl

    wb = openpyxl.load_workbook(path)
    ws = wb[SHEET]
    header = [c.value for c in ws[1]]
    if header[:len(OLD_COLUMNS)] != OLD_COLUMNS:
        sys.exit(f"{SHEET}: unexpected header {header[:len(OLD_COLUMNS)]}")
    template = ws.cell(1, 1)
    for j, name in enumerate(COLUMNS, start=1):
        cell = ws.cell(1, j, name)
        if j > len(OLD_COLUMNS):
            for attr in ("font", "fill", "border", "alignment", "number_format"):
                setattr(cell, attr, copy(getattr(template, attr)))
            ws.column_dimensions[cell.column_letter].width = max(14, len(name) + 2)
    if ws.max_row > 1:
        ws.delete_rows(2, ws.max_row - 1)
    for i, rec in enumerate(rounded(df).itertuples(index=False), start=2):
        for j, v in enumerate(rec, start=1):
            ws.cell(i, j, None if (isinstance(v, float) and v != v) else v)

    if SUMMARY_SHEET in wb.sheetnames:
        del wb[SUMMARY_SHEET]
    ss = wb.create_sheet(SUMMARY_SHEET, index=wb.sheetnames.index(SHEET) + 1)
    for j, name in enumerate(summary.columns, start=1):
        c = ss.cell(1, j, name)
        for attr in ("font", "fill", "border", "alignment"):
            setattr(c, attr, copy(getattr(template, attr)))
    for i, (m, v) in enumerate(summary.itertuples(index=False), start=2):
        ss.cell(i, 1, m)
        ss.cell(i, 2, round(float(v), 4) if isinstance(v, float) else v)
    ss.column_dimensions["A"].width = 58
    ss.column_dimensions["B"].width = 12

    if "README" in wb.sheetnames:  # describe both sheets in the workbook's index
        rd = wb["README"]
        entries = {
            SHEET: ("S8", "Cluster detection and boundary agreement between antiSMASH (regions) and GECCO, per "
                          "genome; computed by bgc_analyses/bgc_consolidation.py from the Reannota_New GFF3s"),
            SUMMARY_SHEET: ("S8 (summary)", "Cohort-level consolidation and overlap metrics derived from S8"),
        }
        by_number = {str(rd.cell(r, 1).value).strip(): r for r in range(2, rd.max_row + 1) if rd.cell(r, 1).value}
        for sheet, (num, text) in entries.items():
            r = by_number.get(num) or rd.max_row + 1
            rd.cell(r, 1, num)
            rd.cell(r, 2, sheet)
            rd.cell(r, 3, text)
    wb.save(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dataset", type=Path, default=DATASET_DIR)
    ap.add_argument("--workbook", type=Path, default=WORKBOOK)
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    ap.add_argument("--no-workbook", action="store_true", help="write the TSVs only")
    a = ap.parse_args(argv)

    df, pairs, summary, problems = build(a.dataset)
    a.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out / "bgc_consolidation_per_genome.tsv", sep="\t", index=False)
    pairs.to_csv(a.out / "bgc_consolidation_pairs.tsv", sep="\t", index=False)
    summary.to_csv(a.out / "bgc_consolidation_summary.tsv", sep="\t", index=False)
    print(f"{len(df)} genomes, {int(df['AntiSMASH Count'].sum())} antiSMASH regions, "
          f"{int(df['GECCO Count'].sum())} GECCO clusters, {len(pairs)} overlapping pairs")
    print(summary.to_string(index=False))
    if problems:
        print(f"\n{len(problems)} check(s) to review:")
        for p in problems:
            print("  " + p)
    else:
        print("\nall checks passed")
    if not a.no_workbook:
        update_workbook(a.workbook, df, summary)
        print(f"updated sheets '{SHEET}' and '{SUMMARY_SHEET}' in {a.workbook}")


if __name__ == "__main__":
    main()
