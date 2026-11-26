import csv
from pathlib import Path

import pytest

from pgxrec.phenotype import c19_phenotype, d6_phenotype
from pgxrec.recommend import recommend

DATA = Path(__file__).resolve().parents[1] / "data" / "processed"
TCAS = {"amitriptyline", "clomipramine", "doxepin", "imipramine", "trimipramine"}   # CYP2C19 x CYP2D6 tables
D6_ONLY = {"nortriptyline", "desipramine", "paroxetine", "fluvoxamine"}               # CYP2D6-only tables


# ---------- CYP2C19 diplotype -> phenotype -------------------------------------------------
@pytest.mark.parametrize("diplotype, expected", [
    ("*1/*1", "CYP2C19 Normal Metabolizer"),
    ("*1/*17", "CYP2C19 Rapid Metabolizer"),
    ("*17/*17", "CYP2C19 Ultrarapid Metabolizer"),
    ("*1/*2", "CYP2C19 Intermediate Metabolizer"),
    ("*2/*2", "CYP2C19 Poor Metabolizer"),
    ("*2/*3", "CYP2C19 Poor Metabolizer"),
    ("*2/*17", "CYP2C19 Intermediate Metabolizer"),   # CPIC: *17 does not compensate for no function
    ("*17/*2", "CYP2C19 Intermediate Metabolizer"),   # allele order must not matter
    ("*1/*9", "CYP2C19 Likely Intermediate Metabolizer"),
    ("*2/*9", "CYP2C19 Likely Poor Metabolizer"),
    ("*1/*12", "CYP2C19 Indeterminate"),              # *12 = uncertain function
    ("*38/*38", "CYP2C19 Normal Metabolizer"),        # *38 = reference-sequence haplotype (A at g.94842866)
    ("*2/*38", "CYP2C19 Intermediate Metabolizer"),
])
def test_c19_phenotype(diplotype, expected):
    assert c19_phenotype(diplotype).phenotype == expected


def test_c19_unknown_allele_is_flagged_not_guessed():
    r = c19_phenotype("*1/*999")
    assert r.phenotype == "CYP2C19 Indeterminate"
    assert r.unknown_alleles == ["*999"]


@pytest.mark.parametrize("bad", ["*1", "1/2", "*1/2", "*1/*2/*3", ""])
def test_bad_diplotype_raises(bad):
    with pytest.raises(ValueError):
        c19_phenotype(bad)


# ---------- CYP2D6 diplotype -> phenotype --------------------------------------------------
@pytest.mark.parametrize("diplotype, phenotype, score", [
    ("*1/*1", "CYP2D6 Normal Metabolizer", "2.0"),
    ("*1/*4", "CYP2D6 Intermediate Metabolizer", "1.0"),   # 2019 update: AS 1.0 = IM
    ("*4/*4", "CYP2D6 Poor Metabolizer", "0.0"),
    ("*2/*41", "CYP2D6 Normal Metabolizer", "1.25"),
    ("*10/*10", "CYP2D6 Intermediate Metabolizer", "0.5"),
    ("*1/*1x2", "CYP2D6 Ultrarapid Metabolizer", "3.0"),
])
def test_d6_phenotype(diplotype, phenotype, score):
    r = d6_phenotype(diplotype)
    assert (r.phenotype, r.activity_score) == (phenotype, score)


def test_d6_allele_order_does_not_matter():
    assert d6_phenotype("*4/*1").phenotype == d6_phenotype("*1/*4").phenotype


def test_d6_unlisted_diplotype_is_indeterminate():
    r = d6_phenotype("*1/*9999")
    assert r.phenotype == "CYP2D6 Indeterminate" and not r.in_table


# ---------- data integrity -----------------------------------------------------------------
def test_every_recommendation_has_valid_action():
    with open(DATA / "recommendations.csv", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert {r["action"] for r in rows} <= {"standard", "dose", "avoid", "none"}
    assert {r["drug"] for r in rows} == TCAS | D6_ONLY | {"citalopram", "escitalopram", "sertraline"}


@pytest.mark.parametrize("c19", ["*1/*1", "*1/*17", "*17/*17", "*1/*2", "*2/*2", "*1/*9", "*2/*9", "*1/*12"])
def test_no_missing_cpic_rows_without_cyp2d6(c19):
    res = recommend(c19)
    assert not [r for r in res.undetermined if r.text.startswith("No CPIC row")]


def test_every_c19_d6_phenotype_pair_has_a_row():
    d6_dips = {"PM": "*4/*4", "IM": "*1/*4", "NM": "*1/*1", "UM": "*1/*1x2", "IND": "*1/*9999"}
    for c19 in ["*1/*1", "*1/*17", "*17/*17", "*1/*2", "*2/*2", "*1/*9", "*2/*9", "*1/*12"]:
        for d6 in d6_dips.values():
            res = recommend(c19, d6)
            assert not [r for r in res.undetermined if r.text.startswith("No CPIC row")], (c19, d6)


# ---------- ranking (CYP2C19 only) ---------------------------------------------------------
def names(group):
    return {r.drug for r in group}


def test_nm_all_drugs_tie_at_zero():
    res = recommend("*1/*1")
    assert len(res.groups) == 1 and res.groups[0][0] == 0
    assert len(res.groups[0][1]) == 8


def test_rm_sertraline_and_citalopram_family_before_tcas():
    g = recommend("*1/*17").groups
    assert [s for s, _ in g] == [0, 2]
    assert names(g[0][1]) == {"sertraline", "citalopram", "escitalopram"}
    assert names(g[1][1]) == TCAS


def test_um_citalopram_family_is_critical_and_last():
    g = recommend("*17/*17").groups
    assert [s for s, _ in g] == [0, 2, 3]
    assert names(g[0][1]) == {"sertraline"}
    assert names(g[2][1]) == {"citalopram", "escitalopram"}


def test_pm_ranking():
    g = recommend("*2/*2").groups
    assert [s for s, _ in g] == [1, 2, 3]
    assert names(g[0][1]) == {"sertraline"}
    assert names(g[1][1]) == TCAS - {"amitriptyline"}
    assert names(g[2][1]) == {"amitriptyline", "citalopram", "escitalopram"}


def test_indeterminate_gives_no_ranking():
    res = recommend("*1/*12")                  # no CYP2D6 either
    assert res.groups == [] and res.notes
    assert len(res.undetermined) == 12         # 8 "No recommendation" + 4 CYP2D6-only drugs that need CYP2D6


def test_tcas_flagged_incomplete_without_cyp2d6_and_sertraline_without_cyp2b6():
    res = recommend("*1/*1")
    by = {r.drug: r for r in res.drugs}
    assert by["amitriptyline"].missing_genes == ["CYP2D6"]
    assert by["sertraline"].missing_genes == ["CYP2B6"]
    assert by["escitalopram"].missing_genes == []


# ---------- ranking with CYP2D6 ------------------------------------------------------------
def score_of(res, drug):
    return {r.drug: r.score for r in res.drugs}[drug]


def test_cyp2d6_extremes_make_amitriptyline_critical_even_when_cyp2c19_normal():
    assert score_of(recommend("*1/*1"), "amitriptyline") == 0
    assert score_of(recommend("*1/*1", "*4/*4"), "amitriptyline") == 3      # D6 PM
    assert score_of(recommend("*1/*1", "*1/*1x2"), "amitriptyline") == 3    # D6 UM
    assert score_of(recommend("*1/*1", "*1/*4"), "amitriptyline") == 1      # D6 IM
    assert score_of(recommend("*1/*1", "*1/*1"), "amitriptyline") == 0      # D6 NM


def test_cyp2d6_does_not_change_citalopram_family():
    for d6 in ["*4/*4", "*1/*4", "*1/*1x2"]:
        assert score_of(recommend("*1/*1", d6), "escitalopram") == score_of(recommend("*1/*1"), "escitalopram")


# ---------- CYP2D6-only drugs (nortriptyline, desipramine, paroxetine, fluvoxamine) -----------
def by_drug(res):
    return {r.drug: r for r in res.drugs + res.undetermined}


def test_d6_only_drugs_are_not_ranked_without_cyp2d6():
    res = recommend("*1/*1")
    assert {r.drug for r in res.drugs} .isdisjoint(D6_ONLY)
    needs = {r.drug for r in res.undetermined if r.text.startswith("Requires CYP2D6")}
    assert needs == D6_ONLY


@pytest.mark.parametrize("d6, expected", [
    ("*4/*4",   {"nortriptyline": 3, "desipramine": 2, "paroxetine": 1, "fluvoxamine": 1}),       # PM
    ("*1/*4",   {"nortriptyline": 1, "desipramine": 1, "paroxetine": 1, "fluvoxamine": 0}),       # IM (AS 1.0)
    ("*1/*1",   {"nortriptyline": 0, "desipramine": 0, "paroxetine": 0, "fluvoxamine": 0}),       # NM
    ("*1/*1x2", {"nortriptyline": 3, "desipramine": 2, "paroxetine": 3, "fluvoxamine": None}),    # UM; fluvoxamine: no evidence
])
def test_d6_only_drug_scores(d6, expected):
    res = recommend("*1/*1", d6)
    got = {d: by_drug(res)[d].score for d in D6_ONLY}
    assert got == expected


def test_d6_only_drugs_do_not_depend_on_cyp2c19():
    ref = {d: by_drug(recommend("*1/*1", "*4/*4"))[d].score for d in D6_ONLY}
    for c19 in ["*1/*17", "*2/*2", "*2/*17", "*1/*12"]:      # includes an Indeterminate CYP2C19 phenotype
        assert {d: by_drug(recommend(c19, "*4/*4"))[d].score for d in D6_ONLY} == ref


def test_nortriptyline_strength_depends_on_activity_score_within_im():
    assert by_drug(recommend("*1/*1", "*1/*4"))["nortriptyline"].strength == "Optional"    # AS 1.0
    assert by_drug(recommend("*1/*1", "*10/*10"))["nortriptyline"].strength == "Moderate"  # AS 0.5


def test_d6_indeterminate_gives_no_recommendation_for_d6_only_drugs():
    res = recommend("*1/*1", "*1/*9999")
    assert all(by_drug(res)[d].score is None for d in D6_ONLY)
