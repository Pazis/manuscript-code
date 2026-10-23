import pytest

import bgc

BAKTA = """##gff-version 3
c1\tBakta\tCDS\t100\t400\t.\t+\t0\tID=X_1;product=hypothetical protein;Dbxref=SO:0001217
c1\tBakta\tCDS\t500\t900\t.\t+\t0\tID=X_2;product=polyketide synthase;Dbxref=SO:0001217,KEGG:K00001,COG:Q,EC:2.3.1.94
c1\tBakta\tCDS\t950\t1200\t.\t-\t0\tID=X_3;product=hypothetical protein
c1\tBakta\tCDS\t5000\t5300\t.\t+\t0\tID=X_4;product=hypothetical protein
##FASTA
>c1
ACGT
"""

ENH = """##gff-version 3
c1\tBiopython\tCDS\t100\t400\t.\t+\t.\tID=X_1;product=Methyltransferase%2C type 11;Dbxref=SO:0001217,InterPro:IPR013216,GO:0008168,COG:COG2226,COG_cat:H
c1\tBiopython\tCDS\t500\t900\t.\t+\t.\tID=X_2;product=polyketide synthase;Dbxref=SO:0001217,KEGG:K00001,COG:Q,COG_cat:Q,EC:EC:2.3.1.94
c1\tBiopython\tCDS\t950\t1200\t.\t-\t.\tID=X_3;product=consensus disorder prediction
c1\tBiopython\tCDS\t5000\t5300\t.\t+\t.\tID=X_4;product=hypothetical protein;Dbxref=COG:S
## Predicted BGC features (antiSMASH / GECCO)
c1\tantiSMASH:7.1.0\tbiosynthetic-gene-cluster\t50\t1000\t.\t.\t.\tID=as:c1_bgc1;product=T1PKS
c1\tantiSMASH:7.1.0\tCDS\t100\t400\t.\t+\t0\tID=as:X_1;Parent=as:c1_bgc1;product=hypothetical protein
c1\tGECCO\tbiosynthetic-gene-cluster\t450\t1300\t.\t.\t.\tID=gecco:c1_cluster_1;product=Polyketide
c1\tGECCO\tCDS\t500\t900\t.\t+\t0\tID=X_2;Parent=gecco:c1_cluster_1
"""


@pytest.fixture
def gffs(tmp_path):
    (tmp_path / "b.gff3").write_text(BAKTA)
    (tmp_path / "e.gff3").write_text(ENH)
    return tmp_path / "b.gff3", tmp_path / "e.gff3"


def test_membership_and_gain(gffs):
    cds, clusters = bgc.cds_table("G", *gffs)
    assert clusters.bgc.tolist() == ["as:c1_bgc1"]  # GECCO cluster ignored
    r = cds.set_index("protein")
    assert r.index.tolist() == ["X_1", "X_2"]  # X_3 ends past the cluster (only in GECCO), X_4 outside
    assert r.loc["X_1", "bgcs"] == "as:c1_bgc1" and r.loc["X_1", "bgc_products"] == "T1PKS"
    assert r.loc["X_1", "recovered"] and r.loc["X_1", "dbxrefs_added"] == 2 and r.loc["X_1", "added_dbs"] == "COG,GO"
    assert not r.loc["X_2", "improved"] and r.loc["X_2", "dbxrefs_lost"] == 0  # COG_cat, COG letter, EC:EC: are not gains
    assert r.loc["X_2", "dbxrefs_bakta"] == r.loc["X_2", "dbxrefs_enhanced"] == 2  # KEGG + EC


def test_mobidb_only_is_not_recovered(gffs, tmp_path):
    (tmp_path / "e2.gff3").write_text(ENH.replace("50\t1000", "50\t1300"))  # now X_3 is inside
    r = bgc.cds_table("G", gffs[0], tmp_path / "e2.gff3")[0].set_index("protein")
    assert r.loc["X_3", "is_hypothetical"] and not r.loc["X_3", "recovered"] and not r.loc["X_3", "improved"]


def test_per_genome_and_db(gffs):
    cds, clusters = bgc.cds_table("G", *gffs)
    pg = bgc.per_genome(cds, clusters.assign(genome="G"), ["G", "EMPTY"]).set_index("genome")
    assert pg.loc["G", ["n_clusters", "n_cds", "n_hypothetical", "n_recovered", "n_improved", "dbxrefs_added"]].tolist() == [1, 2, 1, 1, 1, 2]
    assert pg.loc["G", "pct_hyp_recovered"] == 100 and pg.loc["G", "pct_cds_improved"] == 50
    assert pg.loc["EMPTY", ["n_clusters", "n_cds"]].tolist() == [0, 0]
    db = bgc.per_db(cds).set_index("db")
    assert set(db.index) == {"KEGG", "GO", "COG", "EC"}  # SO, InterPro, COG_cat not counted
    assert db.loc["EC", ["bakta", "enhanced", "added"]].tolist() == [1, 1, 0]
    assert db.loc["GO", "added"] == 1 and db.loc["KEGG", "added"] == 0
    s = bgc.summarise(pg.reset_index()).iloc[0]
    assert s.n_genomes == 2 and s.n_genomes_with_clusters == 1 and s.pct_cds_improved_pooled == 50


def test_missing_cds_fails(gffs, tmp_path):
    (tmp_path / "e2.gff3").write_text(ENH.replace("ID=X_4;", "ID=X_9;"))
    with pytest.raises(ValueError, match="missing from enhanced"):
        bgc.cds_table("G", gffs[0], tmp_path / "e2.gff3")
