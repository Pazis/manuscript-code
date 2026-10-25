"""Tests for bgc_consolidation.py on hand-made GFF3s whose answers can be checked by hand."""
import math

import pytest

import bgc_consolidation as bc

CDS = "contig_1\tBiopython\tCDS\t1\t300\t.\t+\t.\tID=X_00001;product=hypothetical protein\n" \
      "contig_2\tBiopython\tCDS\t1\t300\t.\t+\t.\tID=X_00002;product=hypothetical protein\n"


def gff(tmp_path, rows, name="g.gff3"):
    p = tmp_path / name
    body = "".join(f"{c}\t{src}\tbiosynthetic-gene-cluster\t{s}\t{e}\t.\t.\t.\tID={i};product={prod};tool={src}\n"
                   for c, src, s, e, i, prod in rows)
    p.write_text("##gff-version 3\n" + CDS + "## Predicted BGC features (antiSMASH / GECCO) \n" + body)
    return p


def test_metrics_by_hand(tmp_path):
    # region1 1000-1999 overlaps gecco1 1500-2499 (500 bp) and gecco2 1900-1949 (50 bp, inside gecco1 span)
    # region2 5000-5999 has no GECCO partner; gecco3 sits on contig_2 at the same coordinates (must NOT match)
    p = gff(tmp_path, [
        ("contig_1", "antiSMASH", 1000, 1999, "as:contig_1_region1", "NRPS;T1PKS"),
        ("contig_1", "antiSMASH", 5000, 5999, "as:contig_1_region2", "terpene"),
        ("contig_1", "GECCO", 1500, 2499, "gecco:contig_1_cluster_1", "NRP;Polyketide"),
        ("contig_1", "GECCO", 1900, 1949, "gecco:contig_1_cluster_2", "Unknown"),
        ("contig_2", "GECCO", 5000, 5999, "gecco:contig_2_cluster_1", "NRP"),
    ])
    regions, gecco, contigs = bc.read_clusters(p)
    assert regions[0]["product"] == "NRPS;T1PKS" and gecco[0]["product"] == "NRP;Polyketide"
    row, pairs = bc.genome_metrics("G", regions, gecco)
    assert row["AntiSMASH Count"] == 2 and row["GECCO Count"] == 3
    assert row["Overlapping Clusters"] == 1 and row["GECCO Matched"] == 2
    assert row["Total Overlap (bp)"] == 500           # union of GECCO on contig_1 is 1500-2499
    assert row["AntiSMASH Locus (bp)"] == 2000 and row["GECCO Locus (bp)"] == 2000
    assert math.isclose(row["% AntiSMASH Covered"], 25.0) and math.isclose(row["% GECCO Covered"], 25.0)
    assert math.isclose(row["Inter-tool Consensus"], 0.5)
    assert math.isclose(row["Consolidation Ratio"], 2 / 3)
    assert math.isclose(row["Matched Consolidation Ratio"], 0.5)
    assert row["Regions with >1 GECCO"] == 1 and row["GECCO with >1 Region"] == 0
    assert row["AntiSMASH Only"] == 1 and row["GECCO Only"] == 1
    # Jaccard: 500/1500 and 50/1000
    assert math.isclose(row["Avg Overlap Efficiency"], (500 / 1500 + 50 / 1000) / 2)
    assert row["GECCO Products"].split(", ")[0] == "NRP"
    assert not bc.check_row(row)


def test_no_cross_contig_overlap(tmp_path):
    p = gff(tmp_path, [("contig_1", "antiSMASH", 100, 200, "as:contig_1_region1", "terpene"),
                       ("contig_2", "GECCO", 100, 200, "gecco:contig_2_cluster_1", "Terpene")])
    row, pairs = bc.genome_metrics("G", *bc.read_clusters(p)[:2])
    assert row["Overlapping Clusters"] == 0 and row["Total Overlap (bp)"] == 0 and not pairs


def test_refuses_protocluster_gff(tmp_path):
    p = gff(tmp_path, [("contig_1", "antiSMASH", 100, 200, "as:contig_1_bgc1", "terpene")])
    with pytest.raises(ValueError, match="protocluster"):
        bc.read_clusters(p)


def test_efficiency_never_above_one(tmp_path):
    # identical intervals and nested intervals: Jaccard 1 and < 1, never more
    p = gff(tmp_path, [("contig_1", "antiSMASH", 100, 999, "as:contig_1_region1", "x"),
                       ("contig_1", "GECCO", 100, 999, "gecco:contig_1_cluster_1", "y"),
                       ("contig_1", "GECCO", 200, 300, "gecco:contig_1_cluster_2", "y")])
    row, pairs = bc.genome_metrics("G", *bc.read_clusters(p)[:2])
    assert max(p["jaccard"] for p in pairs) == 1
    assert 0 < row["Avg Overlap Efficiency"] <= 1
    assert row["Total Overlap (bp)"] == 900
