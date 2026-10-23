"""Review candidates for manual curation, and Table 2 from the curated rows."""
import pandas as pd

N_CANDIDATES = 20
MAX_PER_GENOME = 3
MAX_PER_KO = 2  # so one KO (e.g. nifK) cannot dominate the list
CRITERIA = {"hypothetical": 1, "module_change": 2, "ips_confirmed": 1, "requires_review": 2}
ANNOTATION = {"Well-supported", "Plausible-limited", "Conflicting"}
BIOLOGY = {"consistent", "unexpected", "inconsistent"}
CURATION_COLS = ["annotation_consistency", "biological_consistency", "include_in_table2", "notes"]
LIST_COLS = ["interpro", "eggnog_pfams", "modules_affected", "change_category", "criteria_met"]


def score_pairs(pko, proteins, contributions, genomes, ko_names) -> pd.DataFrame:
    """Every new (eggNOG) protein–KO pair with its review score and the criteria it meets."""
    new = pko[pko.in_eggnog].merge(
        proteins[["genome", "protein", "bakta_product", "enhanced_product", "interpro", "eggnog_pfams"]],
        on=["genome", "protein"])
    affected = (contributions.groupby(["genome", "protein", "ko"])
                .agg(modules_affected=("module", lambda m: ",".join(sorted(set(m)))),
                     change_category=("change_category", lambda c: ",".join(sorted(set(c)))),
                     requires_review=("requires_review", "any"))
                .reset_index())
    new = new.merge(affected, on=["genome", "protein", "ko"], how="left")
    new["modules_affected"] = new.modules_affected.fillna("")
    new["change_category"] = new.change_category.fillna("")
    new["requires_review"] = new.requires_review.eq(True)
    met = pd.DataFrame({"hypothetical": new.is_hypothetical, "module_change": new.modules_affected != "",
                        "ips_confirmed": new.ips_confirmed, "requires_review": new.requires_review})
    new["score"] = sum(met[c] * w for c, w in CRITERIA.items())
    new["criteria_met"] = met.apply(lambda r: ",".join(c for c in CRITERIA if r[c]), axis=1)
    new["ko_name"] = new.ko.map(ko_names).fillna("")
    new = new.merge(genomes[["accession", "organism", "group"]], left_on="genome", right_on="accession")
    return new


def select(scored: pd.DataFrame, n=N_CANDIDATES, per_genome=MAX_PER_GENOME, per_ko=MAX_PER_KO) -> pd.DataFrame:
    """Top-scoring pairs, capped per genome and per KO, with every criterion represented at least once."""
    ranked = scored.sort_values(["score", "genome", "protein"], ascending=[False, True, True])
    chosen, per_g, per_k = [], {}, {}

    def take(row):
        if row.Index in chosen or per_g.get(row.genome, 0) >= per_genome or per_k.get(row.ko, 0) >= per_ko:
            return False
        chosen.append(row.Index)
        per_g[row.genome] = per_g.get(row.genome, 0) + 1
        per_k[row.ko] = per_k.get(row.ko, 0) + 1
        return True

    for criterion in CRITERIA:  # seed: best pair for each criterion
        for row in ranked[ranked.criteria_met.str.contains(criterion)].itertuples():
            if take(row):
                break
    for row in ranked.itertuples():
        if len(chosen) >= n:
            break
        take(row)
    return ranked.loc[chosen].sort_values(["score", "genome"], ascending=[False, True])


def candidates_table(selected: pd.DataFrame) -> pd.DataFrame:
    out = selected[["genome", "organism", "group", "protein", "ko", "ko_name", "bakta_product", "enhanced_product",
                    "is_hypothetical", "ips_confirmed", "interpro", "eggnog_pfams", "modules_affected",
                    "change_category", "requires_review", "score", "criteria_met"]].copy()
    # this file is hand-edited in spreadsheets; ";" lists need no quoting, so a sloppy import can't shift columns
    out[LIST_COLS] = out[LIST_COLS].apply(lambda s: s.str.replace(",", ";"))
    out.insert(0, "candidate_id", [f"C{i:02d}" for i in range(1, len(out) + 1)])
    for c in CURATION_COLS:
        out[c] = ""
    return out


def table2(curated: pd.DataFrame) -> pd.DataFrame:
    """Validate the curation and build Table 2 from rows marked include_in_table2 = yes."""
    c = curated.fillna("")
    errors = []
    for col, allowed in [("annotation_consistency", ANNOTATION), ("biological_consistency", BIOLOGY)]:
        bad = c[(c[col] != "") & ~c[col].isin(allowed)]
        errors += [f"{r.candidate_id}: {col}={r[col]!r} not in {sorted(allowed)}" for _, r in bad.iterrows()]
    chosen = c[c.include_in_table2.str.strip().str.lower() == "yes"]
    for _, r in chosen.iterrows():
        for col in ("annotation_consistency", "biological_consistency"):
            if r[col] == "":
                errors.append(f"{r.candidate_id}: selected for Table 2 but {col} is empty")
    if errors:
        raise ValueError("curation problems:\n" + "\n".join(errors))
    source = chosen.ips_confirmed.astype(str).str.lower().map(
        {"true": "eggNOG-mapper; InterProScan-confirmed", "false": "eggNOG-mapper"})
    return pd.DataFrame({
        "genome": chosen.genome, "organism": chosen.organism, "protein": chosen.protein,
        "bakta_product": chosen.bakta_product, "new_ko": chosen.ko, "ko_function": chosen.ko_name,
        "evidence_source": source, "initially_hypothetical": chosen.is_hypothetical,
        "affects_modules": chosen.modules_affected.replace("", "no"),
        "annotation_consistency": chosen.annotation_consistency,
        "biological_consistency": chosen.biological_consistency, "notes": chosen.notes,
    })
