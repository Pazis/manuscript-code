import pandas as pd
import pytest

from config import KEGG_DIR, PANEL
from modules import evaluate, evaluate_panel, kos_in, parse

DEFS = pd.read_csv(KEGG_DIR / "module_definitions.tsv", sep="\t").set_index("module")["definition"].to_dict()
TREES = {m: parse(d) for m, d in DEFS.items()}

NIF = {"K02588", "K02586", "K02591"}
VNF = {"K22896", "K22897", "K22898", "K22899"}
PSII = {"K02703", "K02706", "K02705", "K02704", "K02707", "K02708"}
PSI = {"K02689", "K02690", "K02691", "K02692", "K02693", "K02694"}
CALVIN = {"K00855", "K01601", "K00927", "K00134", "K01803", "K01623", "K03841",
          "K00615", "K01100", "K01807", "K01783"}  # one alternative per block, K01602 omitted


def test_panel_parses_and_is_frozen():
    assert set(TREES) == set(PANEL)


@pytest.mark.parametrize("module", [m for m in PANEL if m != "M00611"])
def test_all_kos_complete_and_none_incomplete(module):
    tree = TREES[module]
    assert evaluate(tree, kos_in(tree))["complete"]
    r = evaluate(tree, set())
    assert r["steps_present"] == 0 and r["blocks_missing"] == r["steps_total"]


def test_nitrogenase_without_optional_k00531():
    assert evaluate(TREES["M00175"], NIF)["complete"]


def test_nitrogenase_alternative_vnf_anf():
    assert evaluate(TREES["M00175"], VNF)["complete"]


def test_nitrogenase_partial_complex_not_complete():
    assert not evaluate(TREES["M00175"], {"K02588", "K02586"})["complete"]


def test_calvin_with_and_without_optional_k01602():
    assert evaluate(TREES["M00165"], CALVIN)["complete"]
    assert evaluate(TREES["M00165"], CALVIN | {"K01602"})["complete"]
    assert not evaluate(TREES["M00165"], (CALVIN - {"K01601"}) | {"K01602"})["complete"]


def test_calvin_block_counting():
    r = evaluate(TREES["M00165"], CALVIN - {"K00927", "K01783"})
    assert (r["steps_total"], r["blocks_missing"], r["status"]) == (11, 2, "2 blocks missing")


def test_glycogen_top_level_or_binds_tighter_than_space():
    # M00854: (K00963 (...)),(K00975 (...)) (K00700,K16149) -> 2 blocks
    assert evaluate(TREES["M00854"], set())["steps_total"] == 2
    assert evaluate(TREES["M00854"], {"K00975", "K00703", "K00700"})["complete"]
    assert not evaluate(TREES["M00854"], {"K00975", "K00700"})["complete"]


def test_composite_m00611():
    assert evaluate_panel(TREES, PSII | CALVIN)["M00611"]["complete"]
    assert evaluate_panel(TREES, PSI | CALVIN)["M00611"]["complete"]
    r = evaluate_panel(TREES, PSII | PSI)["M00611"]
    assert not r["complete"] and r["status"] == "1 block missing"


def test_bad_definition_rejected():
    with pytest.raises(ValueError):
        parse("K00001 (K00002")
