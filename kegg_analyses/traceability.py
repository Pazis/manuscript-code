"""Module change classification, traceability, and the nitrogen, benchmark and TCA tables."""
import pandas as pd

import config
from metrics import state_pairs
from modules import evaluate_panel, kos_in, refs, satisfied

NITROGEN_KOS = {
    "nifH": "K02588", "nifD": "K02586", "nifK": "K02591",
    "vnfD/anfD": "K22896", "vnfK/anfK": "K22897", "vnfG/anfG": "K22898", "vnfH": "K22899",
    "narB": "K00367", "nirA": "K00366",
    "nrtA": "K15576", "nrtB": "K15577", "nrtC": "K15578", "nrtD": "K15579",
}


def category(a: dict, d: dict) -> str:
    if a["complete"] and d["complete"]:
        return "complete in both"
    if d["complete"]:
        return "newly complete"
    if d["steps_present"] > a["steps_present"]:
        return "improved"
    return "unchanged"


def genome_contributions(genome: str, gko: pd.DataFrame, trees: dict):
    """Per module: A/D results and the eggNOG KOs that satisfy blocks newly satisfied in D."""
    kos_a, kos_d = set(state_pairs(gko, "A").ko), set(state_pairs(gko, "D").ko)
    res_a, res_d = evaluate_panel(trees, kos_a), evaluate_panel(trees, kos_d)
    done_a = {m: r["complete"] for m, r in res_a.items()}
    done_d = {m: r["complete"] for m, r in res_d.items()}
    new_kos = kos_d - kos_a
    contrib = {}
    for module in sorted(trees, key=lambda m: bool(refs(trees[m]))):  # referenced modules first
        kos = set()
        for i, block in enumerate(trees[module][1]):
            if satisfied(block, kos_d, done_d) and not satisfied(block, kos_a, done_a):
                kos |= kos_in(block) & new_kos
                for child in refs(block):  # composite: credit the child modules' new KOs
                    kos |= contrib.get(child, set())
        contrib[module] = kos
    return res_a, res_d, contrib


def review_reasons(genome, module, res_d, change, pairs, org_modules) -> list:
    reasons = []
    if change in ("newly complete", "improved") and len(pairs) and not pairs.ips_confirmed.any():
        reasons.append("eggNOG-only, no IPS confirmation")
    expected = config.EXPECTATIONS.get(genome, {}).get(module)
    if expected is not None and expected != res_d["complete"]:
        reasons.append(f"contradicts expectation ({'complete' if expected else 'not complete'} expected)")
    if org_modules is not None and module not in config.COMPOSITE and res_d["complete"] != (module in org_modules):
        reasons.append("disagrees with KEGG organism annotation")
    return reasons


def trace(genomes: pd.DataFrame, pko: pd.DataFrame, proteins: pd.DataFrame, trees: dict, names: dict):
    """Returns (supp_genome_module, module_contributions)."""
    prot = proteins.set_index(["genome", "protein"])
    rows, contrib_rows = [], []
    for g in genomes.itertuples():
        gko = pko[pko.genome == g.accession]
        res_a, res_d, contrib = genome_contributions(g.accession, gko, trees)
        org = config.BENCHMARK_ORGS.get(g.accession)
        org_modules = set((config.KEGG_DIR / f"{org}_modules.txt").read_text().split()) if org else None
        for module in trees:
            a, d = res_a[module], res_d[module]
            change = category(a, d)
            pairs = gko[gko.in_eggnog & gko.ko.isin(contrib[module])]
            reasons = review_reasons(g.accession, module, d, change, pairs, org_modules)
            rows.append({
                "genome": g.accession, "organism": g.organism, "group": g.group,
                "module": module, "module_name": names[module], "composite": module in config.COMPOSITE,
                "baseline_status": a["status"], "baseline_steps": f"{a['steps_present']}/{a['steps_total']}",
                "integrated_status": d["status"], "integrated_steps": f"{d['steps_present']}/{d['steps_total']}",
                "change_category": change, "requires_review": bool(reasons), "review_reasons": "; ".join(reasons),
                "contributing_kos": ",".join(sorted(set(pairs.ko))),
                "contributing_proteins": ",".join(f"{p}:{k}" for p, k in zip(pairs.protein, pairs.ko)),
                "annotation_source": "eggNOG-mapper" if len(pairs) else "",
                "n_ips_confirmed": int(pairs.ips_confirmed.sum()),
                "n_hypothetical_in_bakta": int(pairs.is_hypothetical.sum()),
                "interpretation": "",
            })
            for p in pairs.itertuples():
                contrib_rows.append({"genome": g.accession, "module": module, "ko": p.ko, "protein": p.protein,
                                     "bakta_product": prot.at[(g.accession, p.protein), "bakta_product"],
                                     "is_hypothetical": p.is_hypothetical, "ips_confirmed": p.ips_confirmed,
                                     "change_category": change, "requires_review": bool(reasons)})
    return pd.DataFrame(rows), pd.DataFrame(contrib_rows)


def nitrogen_table(genomes, pko, proteins, supp) -> pd.DataFrame:
    rows = []
    status = supp.set_index(["genome", "module"])
    for g in genomes.itertuples():
        gko, gp = pko[pko.genome == g.accession], proteins[proteins.genome == g.accession]
        row = {"genome": g.accession, "organism": g.organism, "group": g.group}
        for gene, ko in NITROGEN_KOS.items():
            hits = gko[gko.ko == ko]
            row[gene] = ";".join(f"{p}({s})" for p, s in zip(hits.protein, hits.provenance))
        is_n = lambda s: s.str.match(r"(nif|vnf|anf)", case=False)
        row["bakta_nif_genes"] = ",".join(sorted(set(gp.bakta_gene[is_n(gp.bakta_gene.fillna(""))])))
        row["enhanced_nif_names"] = ",".join(sorted(set(gp.enhanced_name[is_n(gp.enhanced_name.fillna(""))])))
        for module in ("M00175", "M00531"):
            s = status.loc[(g.accession, module)]
            row[f"{module}_baseline"] = s.baseline_status
            row[f"{module}_integrated"] = s.integrated_status
            row[f"{module}_new_contributors"] = s.contributing_proteins
        row["interpretation"] = ""
        rows.append(row)
    return pd.DataFrame(rows)


def benchmark(pko, supp, trees) -> pd.DataFrame:
    rows = []
    for genome, org in config.BENCHMARK_ORGS.items():
        ref = set((config.KEGG_DIR / f"{org}_kos.txt").read_text().split())
        listed = set((config.KEGG_DIR / f"{org}_modules.txt").read_text().split())
        ref_res = evaluate_panel(trees, ref)
        gko = pko[pko.genome == genome]
        ka, kd = set(state_pairs(gko, "A").ko), set(state_pairs(gko, "D").ko)
        jac = lambda s: len(s & ref) / len(s | ref)
        for r in supp[supp.genome == genome].itertuples():
            rows.append({"genome": genome, "kegg_org": org, "n_kos_kegg": len(ref),
                         "n_kos_baseline": len(ka), "n_kos_integrated": len(kd),
                         "jaccard_baseline": round(jac(ka), 3), "jaccard_integrated": round(jac(kd), 3),
                         "module": r.module, "status_kegg_kos": ref_res[r.module]["status"],
                         "kegg_lists_module": r.module in listed,
                         "status_baseline": r.baseline_status, "status_integrated": r.integrated_status})
    return pd.DataFrame(rows)


def tca_new(pko, proteins, ko_names) -> pd.DataFrame:
    """New eggNOG KOs in the TCA cycle map (map00020). Qualitative list, no pass/fail."""
    m = pd.read_csv(config.KEGG_DIR / "metabolic_kos.tsv", sep="\t")
    tca = set(m.ko[m["map"] == "map00020"])
    new = pko[pko.in_eggnog & pko.ko.isin(tca)].merge(
        proteins[["genome", "protein", "bakta_product"]], on=["genome", "protein"])
    new.insert(3, "ko_name", new.ko.map(ko_names))
    return new[["genome", "protein", "ko", "ko_name", "bakta_product", "is_hypothetical", "ips_confirmed"]]
