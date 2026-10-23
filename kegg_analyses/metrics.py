"""KO coverage metrics and per-genome module reconstruction.

States: A = Bakta, B = A (InterProScan adds no KOs), C = Bakta ∪ eggNOG, D = C.
validate guarantees A ⊆ D per protein, so gains are D − A.
"""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from modules import evaluate_panel

STATES = {"A": True, "B": True, "C": False, "D": False}  # state -> restrict to Bakta KOs
SETS = {"all": None, "hypothetical": True}


def state_pairs(pko: pd.DataFrame, state: str) -> pd.DataFrame:
    return pko[pko.in_bakta] if STATES[state] else pko


def coverage(proteins: pd.DataFrame, pko: pd.DataFrame, tiers: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for genome, prot in proteins.groupby("genome", sort=False):
        gko = pko[pko.genome == genome]
        for set_name, hyp in SETS.items():
            ids = prot.protein if hyp is None else prot.protein[prot.is_hypothetical]
            sub = gko[gko.protein.isin(ids)]
            for tier, kos in tiers.items():
                q = sub[sub.ko.isin(kos)]
                for state in STATES:
                    s = state_pairs(q, state)
                    rows.append({"genome": genome, "set": set_name, "tier": tier, "state": state,
                                 "N": len(ids), "M": s.protein.nunique(),
                                 "pct": 100 * s.protein.nunique() / len(ids) if len(ids) else np.nan,
                                 "pairs": len(s), "distinct_kos": s.ko.nunique(),
                                 "conflicts": int(s.conflict.sum())})
    long = pd.DataFrame(rows)
    a = long[long.state == "A"].set_index(["genome", "set", "tier"])
    d = long[long.state == "D"].set_index(["genome", "set", "tier"])
    wide = pd.DataFrame({"N": a.N, "M_baseline": a.M, "M_integrated": d.M,
                         "pct_baseline": a.pct, "pct_integrated": d.pct, "gain_pp": d.pct - a.pct,
                         "new_pairs": d.pairs - a.pairs, "new_distinct_kos": d.distinct_kos - a.distinct_kos,
                         "conflicts": d.conflicts}).reset_index()
    return long, wide


def hypothetical_recovery(proteins: pd.DataFrame, tiers: dict) -> pd.DataFrame:
    h = proteins[proteins.is_hypothetical].copy()
    new = h.new_kos.str.split(",").map(lambda ks: set(filter(None, ks)))
    for tier, kos in tiers.items():
        h[f"new_{tier}_ko"] = new.map(lambda s: bool(s & kos))
    h["new_ec"] = [bool(set(filter(None, e.split(","))) - set(filter(None, b.split(","))))
                   for e, b in zip(h.enhanced_ec, h.bakta_ec)]
    h["any_signal"] = h.has_ips | h.has_eggnog
    cols = ["has_ips", "has_eggnog", "any_signal", "named_after", *[f"new_{t}_ko" for t in tiers], "new_ec"]
    out = h.groupby("genome", sort=False)[cols].sum()
    out.insert(0, "n_hypothetical", h.groupby("genome", sort=False).size())
    return out.reset_index()


def summarise(wide: pd.DataFrame, recovery: pd.DataFrame) -> pd.DataFrame:
    """Table 1: median [IQR] across genomes, totals, and paired Wilcoxon baseline vs integrated coverage."""
    rows = []
    for (set_name, tier), g in wide.groupby(["set", "tier"], sort=False):
        diff = g.pct_integrated - g.pct_baseline
        p = wilcoxon(g.pct_integrated, g.pct_baseline).pvalue if (diff != 0).any() else np.nan
        for metric in ["N", "M_baseline", "M_integrated", "pct_baseline", "pct_integrated",
                       "gain_pp", "new_pairs", "new_distinct_kos", "conflicts"]:
            rows.append(stat(f"{metric}", set_name, tier, g[metric],
                             p if metric == "gain_pp" else np.nan))
    for metric in recovery.columns.drop("genome"):
        rows.append(stat(f"hypothetical_{metric}", "hypothetical", "", recovery[metric], np.nan))
    return pd.DataFrame(rows)


def stat(metric, set_name, tier, x, p):
    return {"metric": metric, "set": set_name, "tier": tier, "n_genomes": len(x),
            "median": x.median(), "q1": x.quantile(.25), "q3": x.quantile(.75),
            "total": x.sum() if not metric.startswith("pct") and metric != "gain_pp" else np.nan,
            "wilcoxon_p": p}


def reconstruct(pko: pd.DataFrame, genomes: list, trees: dict) -> pd.DataFrame:
    rows = []
    for genome in genomes:
        gko = pko[pko.genome == genome]
        for state in STATES:
            for module, r in evaluate_panel(trees, set(state_pairs(gko, state).ko)).items():
                rows.append({"genome": genome, "state": state, "module": module, **r})
    return pd.DataFrame(rows)


ALL = ("N", "all", "T2")  # pooled denominators: all proteins, all Bakta-hypotheticals
HYP = ("hypothetical_n_hypothetical", "hypothetical", "")
TABLE1_ROWS = [  # (label, (metric, set, tier), denominator key for a pooled % or None)
    ("Proteins analysed", ALL, None),
    ("Proteins initially hypothetical (Bakta)", HYP, ALL),
    ("Hypotheticals with any new functional evidence (IPS or eggNOG)", ("hypothetical_any_signal", "hypothetical", ""), HYP),
    ("Hypotheticals with InterProScan evidence", ("hypothetical_has_ips", "hypothetical", ""), HYP),
    ("Hypotheticals with eggNOG-mapper evidence", ("hypothetical_has_eggnog", "hypothetical", ""), HYP),
    ("Hypotheticals given an informative product", ("hypothetical_named_after", "hypothetical", ""), HYP),
    ("Hypotheticals with a new metabolic KO (T2)", ("hypothetical_new_T2_ko", "hypothetical", ""), HYP),
    ("Hypotheticals with a new panel-module KO (T1)", ("hypothetical_new_T1_ko", "hypothetical", ""), HYP),
    ("Hypotheticals with a new EC number", ("hypothetical_new_ec", "hypothetical", ""), HYP),
    ("Coverage gain among hypotheticals, percentage points (T2)", ("gain_pp", "hypothetical", "T2"), None),
]
for tier, kind in (("T2", "metabolic"), ("T1", "panel-module")):
    TABLE1_ROWS += [
        (f"Proteins with a {kind} KO, Bakta ({tier})", ("M_baseline", "all", tier), ALL),
        (f"Proteins with a {kind} KO, integrated ({tier})", ("M_integrated", "all", tier), ALL),
        (f"Coverage per genome, Bakta, % ({tier})", ("pct_baseline", "all", tier), None),
        (f"Coverage per genome, integrated, % ({tier})", ("pct_integrated", "all", tier), None),
        (f"Coverage gain, percentage points ({tier})", ("gain_pp", "all", tier), None),
        (f"Additional protein–KO associations ({tier})", ("new_pairs", "all", tier), None),
        (f"Additional distinct KO identifiers per genome ({tier})", ("new_distinct_kos", "all", tier), None),
    ]


def num(x, decimals: int) -> str:
    """Plain number, trailing zeros dropped (5010.0 -> 5010, 838.5 stays 838.5); '' for NaN."""
    if pd.isna(x):
        return ""
    s = f"{x:.{decimals}f}"
    return s.rstrip("0").rstrip(".") if "." in s else s


def table1_manuscript(t1: pd.DataFrame) -> pd.DataFrame:
    """Table 1: median and IQR across genomes, total, and pooled % where it applies, as plain numbers."""
    t = t1.fillna({"tier": ""}).set_index(["metric", "set", "tier"])
    rows = []
    for label, key, denom in TABLE1_ROWS:
        r = t.loc[key]
        d = 2 if key[0].startswith("pct") or key[0] == "gain_pp" else 1
        rows.append({"metric": label,
                     "median": num(r["median"], d), "q1": num(r.q1, d), "q3": num(r.q3, d),
                     # a sum of per-genome distinct KO counts counts a shared KO once per genome, so no total
                     "total": "" if key[0] == "new_distinct_kos" else num(r.total, 0),
                     "pooled_pct": num(100 * r.total / t.loc[denom].total, 1) if denom else "",
                     "wilcoxon_p": "" if pd.isna(r.wilcoxon_p) else f"{r.wilcoxon_p:.2g}"})
    return pd.DataFrame(rows)
