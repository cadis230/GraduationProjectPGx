"""Diplotype -> phenotype translation using the CPIC tables (no ML, pure lookup).

CYP2C19: allele -> function (allele functionality table), then
         (function 1, function 2) -> phenotype (CPIC phenotype table).
         Alleles that are not in the functionality table are treated as
         "Uncertain function" and reported in `unknown_alleles`.
CYP2D6:  the diplotype is looked up directly in CPIC's pre-enumerated
         diplotype table (~16.8k rows, includes copy-number notation such as
         *1/*2x2). Diplotypes not in the table are "Indeterminate".
"""
import csv
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data" / "processed"

UNCERTAIN = "Uncertain function"
C19_INDETERMINATE = "CYP2C19 Indeterminate"
D6_INDETERMINATE = "CYP2D6 Indeterminate"


@dataclass
class PhenotypeResult:
    gene: str
    diplotype: str
    phenotype: str
    activity_score: str = ""            # CYP2D6 only
    functions: tuple = ()               # CYP2C19 only: function of each allele
    unknown_alleles: list = field(default_factory=list)
    in_table: bool = True

    @property
    def indeterminate(self) -> bool:
        return self.phenotype.endswith("Indeterminate")


def split_diplotype(diplotype: str) -> tuple:
    parts = diplotype.replace(" ", "").split("/")
    if len(parts) != 2 or not all(p.startswith("*") and len(p) > 1 for p in parts):
        raise ValueError(f"Diplotype must look like '*1/*2', got {diplotype!r}")
    return parts[0], parts[1]


@lru_cache(maxsize=1)
def _c19_tables():
    with open(DATA / "c19_allele_function.csv", newline="", encoding="utf-8") as fh:
        alleles = {r["allele"]: r["function"] for r in csv.DictReader(fh)}
    with open(DATA / "c19_phenotype_rules.csv", newline="", encoding="utf-8") as fh:
        rules = {(r["function_1"], r["function_2"]): r["phenotype"] for r in csv.DictReader(fh)}
    return alleles, rules


@lru_cache(maxsize=1)
def _d6_table():
    with open(DATA / "d6_diplotype_phenotype.csv", newline="", encoding="utf-8") as fh:
        return {r["diplotype"]: (r["phenotype"], r["activity_score"]) for r in csv.DictReader(fh)}


def c19_phenotype(diplotype: str) -> PhenotypeResult:
    alleles, rules = _c19_tables()
    a1, a2 = split_diplotype(diplotype)
    unknown = [a for a in (a1, a2) if a not in alleles]
    f1, f2 = (alleles.get(a, UNCERTAIN) for a in (a1, a2))
    phenotype = rules.get((f1, f2), C19_INDETERMINATE)
    return PhenotypeResult("CYP2C19", f"{a1}/{a2}", phenotype,
                           functions=(f1, f2), unknown_alleles=unknown,
                           in_table=(f1, f2) in rules)


def d6_phenotype(diplotype: str) -> PhenotypeResult:
    table = _d6_table()
    a1, a2 = split_diplotype(diplotype)
    for key in (f"{a1}/{a2}", f"{a2}/{a1}"):   # table stores one allele order
        if key in table:
            phenotype, score = table[key]
            return PhenotypeResult("CYP2D6", key, phenotype, activity_score=score)
    return PhenotypeResult("CYP2D6", f"{a1}/{a2}", D6_INDETERMINATE, in_table=False)
