import pandas as pd
import pytest

import review as rv
import run_kegg_eval
from test_modules import CALVIN, PSII, TREES
from traceability import category, genome_contributions


def profile(bakta: set, eggnog: set) -> pd.DataFrame:
    rows = [(k, True) for k in bakta] + [(k, False) for k in eggnog]
    return pd.DataFrame(rows, columns=["ko", "in_bakta"])


def test_new_ko_credited_to_the_block_it_completes_and_to_the_composite():
    gko = profile(PSII | CALVIN - {"K00927"}, {"K00927", "K99999"})
    res_a, res_d, contrib = genome_contributions("g", gko, TREES)
    assert contrib["M00165"] == {"K00927"}  # K99999 is in no block
    assert contrib["M00611"] == {"K00927"}  # via its child M00165
    assert category(res_a["M00165"], res_d["M00165"]) == "newly complete"
    assert contrib["M00161"] == set()


def test_category_improved_and_unchanged():
    a = {"complete": False, "steps_present": 3}
    assert category(a, {"complete": False, "steps_present": 4}) == "improved"
    assert category(a, dict(a)) == "unchanged"


def scored(n_genomes=10, per=5) -> pd.DataFrame:
    rows = []
    for g in range(n_genomes):
        for p in range(per):
            crit = ["module_change"] + (["requires_review"] if (g, p) == (9, 4) else [])
            rows.append({"genome": f"G{g}", "protein": f"P{g}_{p}", "ko": f"K{g:02d}{p:03d}",
                         "score": len(crit), "criteria_met": ",".join(crit)})
    rows[0]["ko"] = rows[1]["ko"] = rows[2]["ko"] = "K00001"
    return pd.DataFrame(rows)


def test_selection_caps_and_covers_every_present_criterion():
    sel = rv.select(scored())
    assert len(sel) == rv.N_CANDIDATES
    assert sel.groupby("genome").size().max() <= rv.MAX_PER_GENOME
    assert sel.groupby("ko").size().max() <= rv.MAX_PER_KO
    assert sel.criteria_met.str.contains("requires_review").any()  # rare criterion still seeded


def test_candidates_table_lists_use_semicolons():
    sel = pd.DataFrame([{"genome": "G", "organism": "O", "group": "X", "protein": "P", "ko": "K00001",
                         "ko_name": "a, b; c", "bakta_product": "", "enhanced_product": "", "is_hypothetical": True,
                         "ips_confirmed": True, "interpro": "IPR1,IPR2", "eggnog_pfams": "A,B",
                         "modules_affected": "M1,M2", "change_category": "improved,newly complete",
                         "requires_review": False, "score": 4, "criteria_met": "hypothetical,module_change"}])
    out = rv.candidates_table(sel)
    assert not out[rv.LIST_COLS].apply(lambda s: s.str.contains(",")).any().any()
    assert out.criteria_met.iloc[0] == "hypothetical;module_change"


def curated(**overrides) -> pd.DataFrame:
    row = {"candidate_id": "C01", "genome": "G", "organism": "O", "protein": "P", "ko": "K00001", "ko_name": "x",
           "bakta_product": "hypothetical protein", "is_hypothetical": "True", "ips_confirmed": "True",
           "modules_affected": "M00165", "annotation_consistency": "Well-supported",
           "biological_consistency": "consistent", "include_in_table2": "yes", "notes": ""}
    row.update(overrides)
    return pd.DataFrame([row, {**row, "candidate_id": "C02", "include_in_table2": ""}])


def test_table2_keeps_selected_rows():
    t = rv.table2(curated())
    assert len(t) == 1 and t.evidence_source.iloc[0] == "eggNOG-mapper; InterProScan-confirmed"


@pytest.mark.parametrize("override", [{"annotation_consistency": ""}, {"biological_consistency": "fine"}])
def test_table2_rejects_incomplete_or_misspelled_curation(override):
    with pytest.raises(ValueError):
        rv.table2(curated(**override))


def test_review_never_overwrites_curation(tmp_path, monkeypatch):
    target = tmp_path / "review_candidates.csv"
    target.write_text("curated by hand\n")
    monkeypatch.setattr(run_kegg_eval, "CANDIDATES", target)
    monkeypatch.setattr(run_kegg_eval, "OUT", tmp_path)
    monkeypatch.setattr(run_kegg_eval, "load", lambda name: pd.DataFrame())
    monkeypatch.setattr(rv, "score_pairs", lambda *a: scored())
    monkeypatch.setattr(rv, "candidates_table", lambda sel: sel)
    run_kegg_eval.review()
    assert target.read_text() == "curated by hand\n"
    assert (tmp_path / "review_candidates_new.csv").exists()


def test_table1_numbers_are_plain_and_unrounded():
    from metrics import num
    assert [num(40, 0), num(5010.0, 1), num(838.5, 1), num(27.0, 1), num(float("nan"), 1)] == ["40", "5010", "838.5", "27", ""]
