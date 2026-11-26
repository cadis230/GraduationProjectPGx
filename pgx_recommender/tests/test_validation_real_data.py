"""Regression test on REAL data: GeT-RM consensus genotypes vs NYGC 1000 Genomes GRCh38 phased VCF (CYP2C19 region)."""
import pytest
from pgxrec.validate import run


@pytest.fixture(scope="module")
def result():
    return run()


def test_vcf_matches_reference_build(result):
    assert result[1]["ref_mismatch_positions"] == []       # REF bases agree with the CPIC table (GRCh38)


def test_all_samples_phased(result):
    assert result[1]["unphased"] == 0


def test_concordance_on_real_samples(result):
    d, s = result
    assert s["n_in_vcf"] == 88
    assert s["diplotype_ok"] == s["phenotype_ok"] == s["ranking_ok"] == 88
