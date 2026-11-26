"""Phenotype -> per-drug CPIC recommendation -> ordinal score -> ranking.

Each drug gets exactly ONE CPIC recommendation per patient (CPIC tables already
combine genes in one row), so each drug gets exactly one score:

    0  standard dose
    1  dose adjustment
    2  avoid / switch, recommendation strength "Optional"
    3  avoid / switch, strength "Moderate" or "Strong"   (= "critical")
    None  CPIC gives no recommendation (reported separately, never ranked)

Lower is better. Ties are reported as ties, not broken artificially.
The score measures genetic *compatibility risk* according to CPIC, not efficacy.
"""
import csv
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Optional

from .phenotype import PhenotypeResult, c19_phenotype, d6_phenotype

DATA = Path(__file__).resolve().parents[1] / "data" / "processed"

SCORE_LABELS = {
    0: "standard dose",
    1: "dose adjustment",
    2: "avoid / switch (optional)",
    3: "avoid / switch (moderate-strong) - CRITICAL",
}
D6_NO_RESULT = "CYP2D6 No Result"
B6_NO_RESULT = "CYP2B6 No Result"   # CYP2B6 is not supported yet
REQUIRES_D6 = "Requires CYP2D6 genotype (the CPIC table has no No-Result fallback)"


def score(action: str, strength: str) -> Optional[int]:
    if action == "none":
        return None
    if action == "standard":
        return 0
    if action == "dose":
        return 1
    if action == "avoid":
        return 3 if strength in ("Moderate", "Strong") else 2
    raise ValueError(f"unknown action {action!r}")


@dataclass
class DrugRecommendation:
    drug: str
    drug_class: str
    score: Optional[int]
    action: str
    strength: str
    text: str
    missing_genes: list = field(default_factory=list)  # genes the CPIC table could use but we did not supply


@dataclass
class Result:
    cyp2c19: PhenotypeResult
    cyp2d6: Optional[PhenotypeResult]
    drugs: list                        # DrugRecommendation, sorted by score
    groups: list                       # [(score, [DrugRecommendation, ...])] best first
    undetermined: list                 # DrugRecommendation with score None / no CPIC row
    notes: list = field(default_factory=list)


def _as_key(d6) -> str:
    """Canonical CYP2D6 activity score string, or 'ANY' when it is not a plain number (e.g. '>=3.0')."""
    try:
        return str(float(d6.activity_score))
    except (TypeError, ValueError):
        return "ANY"


@lru_cache(maxsize=1)
def _load():
    with open(DATA / "recommendations.csv", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    table = {(r["drug"], r["cyp2c19"], r["cyp2d6"], r["cyp2b6"], r["cyp2d6_as"]): r for r in rows}
    drugs = {}
    for r in rows:
        d = drugs.setdefault(r["drug"], {"class": r["drug_class"], "c19": False, "d6": False,
                                         "b6": False, "d6_fallback": False})
        d["c19"] |= r["cyp2c19"] != "ANY"
        d["d6"] |= r["cyp2d6"] != "ANY"
        d["b6"] |= r["cyp2b6"] != "ANY"
        d["d6_fallback"] |= r["cyp2d6"] == D6_NO_RESULT
    return table, drugs


def recommend(cyp2c19_diplotype: str, cyp2d6_diplotype: Optional[str] = None) -> Result:
    table, drugs = _load()
    c19 = c19_phenotype(cyp2c19_diplotype)
    d6 = d6_phenotype(cyp2d6_diplotype) if cyp2d6_diplotype else None
    notes = []
    if c19.unknown_alleles:
        notes.append(f"CYP2C19 alleles not in CPIC functionality table (treated as uncertain function): "
                     f"{', '.join(c19.unknown_alleles)}")
    if c19.indeterminate:
        notes.append("CYP2C19 phenotype is Indeterminate: CPIC gives no recommendation "
                     "(this is where predicted allele function from the novel-variant module would plug in).")
    if d6 is not None and d6.indeterminate:
        notes.append("CYP2D6 phenotype is Indeterminate; CYP2D6-dependent recommendations are limited.")

    recs = []
    for drug, meta in sorted(drugs.items()):
        if meta["d6"] and d6 is None and not meta["d6_fallback"]:
            recs.append(DrugRecommendation(drug, meta["class"], None, "none", "", REQUIRES_D6, ["CYP2D6"]))
            continue
        c19_key = c19.phenotype if meta["c19"] else "ANY"
        d6_key = (d6.phenotype if d6 else D6_NO_RESULT) if meta["d6"] else "ANY"
        b6_key = B6_NO_RESULT if meta["b6"] else "ANY"
        as_key = _as_key(d6) if (meta["d6"] and d6) else "ANY"
        missing = []
        if meta["d6"] and d6 is None:
            missing.append("CYP2D6")
        if meta["b6"]:
            missing.append("CYP2B6")
        row = (table.get((drug, c19_key, d6_key, b6_key, as_key))
               or table.get((drug, c19_key, d6_key, b6_key, "ANY")))
        if row is None:
            recs.append(DrugRecommendation(drug, meta["class"], None, "none", "", "No CPIC row for this combination", missing))
            continue
        recs.append(DrugRecommendation(drug, meta["class"], score(row["action"], row["strength"]),
                                       row["action"], row["strength"], row["text"], missing))

    ranked = sorted((r for r in recs if r.score is not None), key=lambda r: (r.score, r.drug))
    undetermined = [r for r in recs if r.score is None]
    groups = []
    for r in ranked:
        if groups and groups[-1][0] == r.score:
            groups[-1][1].append(r)
        else:
            groups.append((r.score, [r]))
    return Result(c19, d6, ranked, groups, undetermined, notes)


def format_result(res: Result) -> str:
    lines = [f"CYP2C19 {res.cyp2c19.diplotype}: {res.cyp2c19.phenotype}"]
    if res.cyp2d6:
        lines.append(f"CYP2D6  {res.cyp2d6.diplotype}: {res.cyp2d6.phenotype} (activity score {res.cyp2d6.activity_score})")
    else:
        lines.append("CYP2D6  not provided -> CPIC 'No Result' rows used for TCAs")
    for n in res.notes:
        lines.append(f"NOTE: {n}")
    lines.append("")
    for i, (sc, grp) in enumerate(res.groups, 1):
        lines.append(f"{i}. score {sc} - {SCORE_LABELS[sc]}")
        for r in grp:
            flag = f"  [incomplete: no {'/'.join(r.missing_genes)}]" if r.missing_genes else ""
            lines.append(f"     {r.drug:<14}{r.drug_class:<5} strength={r.strength}{flag}")
    if res.undetermined:
        why = lambda r: "needs CYP2D6" if r.text == REQUIRES_D6 else "no CPIC recommendation"
        lines.append("Not ranked: " + ", ".join(f"{r.drug} ({why(r)})" for r in res.undetermined))
    return "\n".join(lines)
