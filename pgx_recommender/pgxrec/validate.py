"""Validate VCF -> diplotype -> phenotype -> ranking against GeT-RM consensus genotypes.
    python -m pgxrec.validate
"""
import json
from pathlib import Path

import pandas as pd

from .matcher import call_diplotype, load_definitions, numkey, read_vcf
from .phenotype import c19_phenotype
from .recommend import recommend

V = Path(__file__).resolve().parents[1] / "data" / "validation"


def norm(dip):            # GeT-RM does not distinguish *38 (reference) from *1
    return "/".join(sorted((a.replace("*38", "*1") for a in dip.split("/")), key=numkey))


def groups(dip):
    return [(s, sorted(r.drug for r in g)) for s, g in recommend(dip).groups]


def run(vcf=V / "cyp2c19_getrm.vcf.gz", truth=V / "getrm_cyp2c19_truth.csv"):
    alleles, refbase = load_definitions()
    samples, haps, mism = read_vcf(vcf, refbase)
    t = pd.read_csv(truth)
    all_conf = t[t.confirmed_standard & t.in_nygc_3202]
    t = all_conf[all_conf.sample_id.isin(samples)]
    rows = []
    for _, r in t.iterrows():
        h1, h2, flags = haps[r.sample_id]
        dip, tied = call_diplotype(h1, h2, alleles, refbase)
        ok_dip = "?" not in dip
        pt, tt = (c19_phenotype(norm(dip)) if ok_dip else None), c19_phenotype(r.cyp2c19_diplotype)
        rows.append(dict(sample=r.sample_id, truth=r.cyp2c19_diplotype, predicted=dip,
                         truth_phenotype=tt.phenotype.replace("CYP2C19 ", ""), pred_phenotype=pt.phenotype.replace("CYP2C19 ", "") if pt else "Unresolved",
                         diplotype_ok=norm(dip) == r.cyp2c19_diplotype, phenotype_ok=bool(pt) and pt.phenotype == tt.phenotype,
                         ranking_ok=ok_dip and groups(norm(dip)) == groups(r.cyp2c19_diplotype), tied=tied, flags=",".join(sorted(flags))))
    d = pd.DataFrame(rows)
    d.to_csv(V / "results.csv", index=False)
    s = dict(n_truth_in_nygc=int(len(all_conf)), n_in_vcf=int(len(t)), n_missing_from_vcf=int(len(all_conf) - len(t)),
             ref_mismatch_positions=mism, diplotype_ok=int(d.diplotype_ok.sum()), phenotype_ok=int(d.phenotype_ok.sum()),
             ranking_ok=int(d.ranking_ok.sum()), tied_calls=int(d.tied.sum()), unphased=int(d["flags"].str.contains("unphased").sum()),
             by_phenotype=d.groupby("truth_phenotype").agg(n=("phenotype_ok", "size"), phenotype_ok=("phenotype_ok", "sum")).to_dict("index"))
    (V / "summary.json").write_text(json.dumps(s, indent=1, default=int), encoding="utf-8")
    return d, s


if __name__ == "__main__":
    d, s = run()
    print(json.dumps(s, indent=1, default=int))
    bad = d[~d.diplotype_ok]
    print("\nUyuşmayan diplotipler:", len(bad))
    print(bad[["sample", "truth", "predicted", "truth_phenotype", "pred_phenotype", "phenotype_ok", "tied"]].to_string(index=False))
