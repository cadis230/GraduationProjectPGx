"""Convert the raw CPIC / ClinPGx Excel tables in data/raw into clean CSVs.

Run once (or whenever the raw tables are updated):

    python -m pgxrec.build_data

Outputs (data/processed):
    c19_allele_function.csv     allele -> clinical function
    c19_phenotype_rules.csv     (function 1, function 2) -> phenotype  (both orders)
    d6_diplotype_phenotype.csv  diplotype -> activity score, phenotype
    recommendations.csv         one row per (drug, CYP2C19, CYP2D6, CYP2B6 phenotype)

The recommendation text -> action mapping lives in config/action_map.csv and is
meant to be reviewed by a human. If a text is not in the map the build fails,
so a guideline update can never be silently misclassified.
"""
import csv
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
ACTION_MAP = ROOT / "config" / "action_map.csv"

DRUGS = {
    "amitriptyline": "TCA",
    "clomipramine": "TCA",
    "doxepin": "TCA",
    "imipramine": "TCA",
    "trimipramine": "TCA",
    "citalopram": "SSRI",
    "escitalopram": "SSRI",
    "sertraline": "SSRI",
    "nortriptyline": "TCA",
    "desipramine": "TCA",
    "paroxetine": "SSRI",
    "fluvoxamine": "SSRI",
}
ACTIONS = {"standard", "dose", "avoid", "none"}


def norm(x) -> str:
    return re.sub(r"\s+", " ", str(x).strip())


def norm_as(x) -> str:
    """Activity score as a canonical string ('1.0', '0.75'); non-numeric values (e.g. '>=3.0') stay as text."""
    try:
        return str(float(x))
    except (TypeError, ValueError):
        return norm(x)


def is_numeric(x: str) -> bool:
    try:
        float(x)
        return True
    except ValueError:
        return False


def propose_action(text: str) -> str:
    """Rule-based *proposal* used only to pre-fill config/action_map.csv."""
    t = norm(text).lower()
    if t.startswith("no recommendation"):
        return "none"
    if (
        t.startswith("avoid")
        or t.startswith("consider alternative drug")
        or re.match(r"select (an )?alternative", t)
        or re.match(r"consider a clinically appropriate (alternative )?antidepressant not", t)
    ):
        return "avoid"
    if t.startswith("initiate therapy with recommended starting dose. if patient"):
        return "standard"
    if (
        t.startswith("consider a lower starting dose")
        or re.match(r"consider a \d+% reduction", t)
        or "slower titration" in t
    ):
        return "dose"
    if t.startswith("initiate therapy with recommended starting dose"):
        return "standard"
    raise ValueError(f"Cannot classify recommendation text: {text!r}")


def build_c19_alleles() -> None:
    df = pd.read_excel(RAW / "CYP2C19_allele_functionality_reference.xlsx",
                       sheet_name="Allele Function", header=1)
    assert "Allele Clinical Functional Status" in df.columns[3], df.columns[3]
    out = pd.DataFrame({"allele": df.iloc[:, 0].map(norm), "function": df.iloc[:, 3].map(norm)})
    out = out[out["allele"].str.startswith("*")].drop_duplicates()
    assert out["allele"].is_unique
    out.to_csv(OUT / "c19_allele_function.csv", index=False)
    print(f"c19_allele_function.csv: {len(out)} alleles")


def build_c19_rules() -> None:
    df = pd.read_excel(RAW / "CYP2C19_phenotypes.xlsx", sheet_name="Phenotype", header=1)
    rows = set()
    for _, r in df.iterrows():
        f1, f2, ph = norm(r["Allele 1 Function"]), norm(r["Allele 2 Function"]), norm(r["Phenotype"])
        rows.add((f1, f2, ph))
        rows.add((f2, f1, ph))
    out = pd.DataFrame(sorted(rows), columns=["function_1", "function_2", "phenotype"])
    clash = out.groupby(["function_1", "function_2"])["phenotype"].nunique()
    assert (clash == 1).all(), "conflicting phenotype rules"
    out.to_csv(OUT / "c19_phenotype_rules.csv", index=False)
    print(f"c19_phenotype_rules.csv: {len(out)} rows")


def build_d6_diplotypes() -> None:
    df = pd.read_excel(RAW / "CYP2D6_Diplotype_Phenotype_Table.xlsx", sheet_name="Diplotypes")
    out = pd.DataFrame({
        "diplotype": df.iloc[:, 0].map(norm),
        "activity_score": df.iloc[:, 1].map(norm),
        "phenotype": df.iloc[:, 2].map(norm),
    })
    assert out["diplotype"].is_unique
    out.to_csv(OUT / "d6_diplotype_phenotype.csv", index=False)
    print(f"d6_diplotype_phenotype.csv: {len(out)} diplotypes")


def load_action_map(texts: set) -> dict:
    if not ACTION_MAP.exists():
        with open(ACTION_MAP, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["text", "action"])
            for t in sorted(texts):
                w.writerow([t, propose_action(t)])
        print(f"Wrote proposed {ACTION_MAP.name} - REVIEW IT before trusting results.")
    with open(ACTION_MAP, newline="", encoding="utf-8") as fh:
        amap = {r["text"]: r["action"] for r in csv.DictReader(fh)}
    bad = {a for a in amap.values() if a not in ACTIONS}
    assert not bad, f"unknown actions in action_map: {bad}"
    missing = texts - set(amap)
    if missing:
        raise SystemExit(f"{len(missing)} recommendation texts missing from action_map.csv, e.g. {sorted(missing)[0]!r}")
    return amap


def build_recommendations() -> None:
    frames = []
    for drug, drug_class in DRUGS.items():
        df = pd.read_excel(RAW / f"{drug}_recommendation.xlsx", sheet_name=0)
        df.columns = [norm(c) for c in df.columns]
        out = pd.DataFrame({
            "drug": drug,
            "drug_class": drug_class,
            "cyp2c19": df["CYP2C19 Phenotype"].map(norm) if "CYP2C19 Phenotype" in df else "ANY",
            "cyp2d6": df["CYP2D6 Phenotype"].map(norm) if "CYP2D6 Phenotype" in df else "ANY",
            "cyp2b6": df["CYP2B6 Phenotype"].map(norm) if "CYP2B6 Phenotype" in df else "ANY",
            "cyp2d6_as": df["CYP2D6 Activity Score"].map(norm_as) if "CYP2D6 Activity Score" in df else "ANY",
            "strength": df["Classification of Recommendation"].map(norm),
            "text": df["Therapeutic Recommendation"].map(norm),
        })
        # Activity-score rows normally collapse onto one row per phenotype combination.
        # If a phenotype combination has conflicting recommendations across activity scores
        # (happens for nortriptyline CYP2D6 IM), keep the activity score in the key instead.
        key = ["drug", "cyp2c19", "cyp2d6", "cyp2b6"]
        parts = []
        for _, g in out.groupby(key, sort=False):
            if g[["strength", "text"]].drop_duplicates().shape[0] == 1:
                parts.append(g.assign(cyp2d6_as="ANY").drop_duplicates())
            else:
                assert g["cyp2d6_as"].map(is_numeric).all(), f"{drug}: conflicting rows without numeric activity score\n{g}"
                parts.append(g.drop_duplicates())
        out = pd.concat(parts, ignore_index=True)
        dup = out[out.duplicated(key + ["cyp2d6_as"], keep=False)]
        assert dup.empty, f"{drug}: conflicting recommendations for one phenotype combination\n{dup}"
        frames.append(out)
    rec = pd.concat(frames, ignore_index=True)
    amap = load_action_map(set(rec["text"]))
    rec["action"] = rec["text"].map(amap)
    rec.to_csv(OUT / "recommendations.csv", index=False)
    print(f"recommendations.csv: {len(rec)} rows, {rec['drug'].nunique()} drugs, "
          f"{rec['text'].nunique()} distinct texts")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_c19_alleles()
    build_c19_rules()
    build_d6_diplotypes()
    build_recommendations()


if __name__ == "__main__":
    main()
