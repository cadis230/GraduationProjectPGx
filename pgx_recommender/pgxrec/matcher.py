"""GRCh38 phased VCF -> CYP2C19 diplotype, using the CPIC allele definition table (no ML).

Per haplotype: every allele whose defined positions are all satisfied competes (IUPAC cells accept either
base); the allele defining the most positions wins. *38 = reference sequence, *1 = G at g.94842866 only.
Positions absent from the VCF are hom-ref (joint-called cohort VCF only lists variable sites).
Ambiguous IUPAC cells in the table (R/Y/M) are treated as "either base", i.e. not defining.
"""
import gzip
import re
from functools import lru_cache
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFINITION = ROOT / "data" / "raw" / "CYP2C19_allele_definition_table.xlsx"
DEFAULT, REFERENCE = "*1", "*38"


def numkey(a: str):
    m = re.match(r"\*(\d+)(.*)", a)
    return (int(m.group(1)), m.group(2))


IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT"}


@lru_cache(maxsize=1)
def load_definitions():
    """-> ({allele: {pos: allowed_bases}}, {pos: reference_base})  (GRCh38, NC_000010.11)
    Every non-blank cell of the CPIC table is a defined position; IUPAC codes (R, Y, M...) allow either base."""
    raw = pd.read_excel(DEFINITION, sheet_name=0, header=None)
    cols = {}
    for j in range(1, raw.shape[1]):
        m = re.match(r"g\.(\d+)([ACGT])>([ACGT])$", str(raw.iloc[3, j]))
        if m:
            cols[j] = (int(m.group(1)), m.group(2))
    alleles = {}
    for i in range(7, raw.shape[0]):
        name = str(raw.iloc[i, 0]).strip()
        if name.startswith("*"):
            alleles[name] = {pos: set(IUPAC[str(raw.iloc[i, j]).strip()]) for j, (pos, ref) in cols.items()
                             if pd.notna(raw.iloc[i, j]) and str(raw.iloc[i, j]).strip() in IUPAC}
    return alleles, {pos: ref for pos, ref in cols.values()}


def read_vcf(path, refbase):
    """-> (samples, {sample: (hap1, hap2, flags)}, reference_mismatch_positions); haps map pos -> base"""
    opener = gzip.open if str(path).endswith(".gz") else open
    haps, samples, mism = {}, [], set()
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = p[9:]
                haps = {s: ({}, {}, set()) for s in samples}
                continue
            pos, ref, alts = int(p[1]), p[3], p[4].split(",")
            if pos not in refbase or len(ref) != 1 or any(len(a) != 1 for a in alts):
                continue
            if ref != refbase[pos]:
                mism.add(pos)
                continue
            bases = [ref] + alts
            for s, cell in zip(samples, p[9:]):
                gt = cell.split(":")[0]
                sep = "|" if "|" in gt else "/"
                for k, idx in enumerate(gt.split(sep)[:2]):
                    if idx == ".":
                        haps[s][2].add("missing_gt")
                        continue
                    b = bases[int(idx)]
                    if b != ref or pos not in haps[s][k]:
                        haps[s][k][pos] = b
                if sep == "/":
                    haps[s][2].add("unphased")
    return samples, haps, sorted(mism)


FUNCTION_RANK = {"No function": 0, "Decreased function": 1, "Uncertain function": 2, "Normal function": 3, "Increased function": 4}


def call_haplotype(hap, alleles, refbase):
    """-> (allele, tied_alleles or None, matched)
    All alleles whose defined positions are satisfied compete; the one defining the most positions wins.
    Exact ties are flagged and resolved conservatively (lowest enzyme function first, then lowest number)."""
    from .phenotype import _c19_tables
    funcs = _c19_tables()[0]
    hits = [(len(cells), n) for n, cells in alleles.items()
            if all(hap.get(p, refbase[p]) in ok for p, ok in cells.items())]
    if not hits:
        return "?", None
    top = max(h[0] for h in hits)
    best = sorted((n for l, n in hits if l == top),
                  key=lambda n: (FUNCTION_RANK.get(funcs.get(n, "Uncertain function"), 2), numkey(n)))
    return best[0], (sorted(best, key=numkey) if len(best) > 1 else None)


def call_diplotype(h1, h2, alleles, refbase):
    (a1, t1), (a2, t2) = call_haplotype(h1, alleles, refbase)[:2], call_haplotype(h2, alleles, refbase)[:2]
    return "/".join(sorted((a1, a2), key=lambda a: numkey(a) if a != "?" else (999, ""))), bool(t1 or t2)
