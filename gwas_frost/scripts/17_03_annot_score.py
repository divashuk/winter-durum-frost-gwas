# -*- coding: utf-8 -*-
"""Раздел 17, часть A, шаг 3. Аннотация (Ensembl Plants BioMart) и баллы генов новых интервалов I09–I12 —
scripts/14_biomart.py и scripts/14_score_genes.py без изменения метода (исполняются из исходного текста с заменой путей:
genes14 -> genes17; список интервалов -> tables/17_gene_intervals.csv; выходные таблицы -> WSL genes17/score_*.csv).
Затем — объединение с генами интервалов раздела 14 (7B, 1B, 2B 767, 2B 776, 4A: tables/14_genes_all_intervals.csv) -> WSL genes17/genes_all9.csv."""
import os, sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
G = r"/work/scratch/genes17"


def run(script, subs):
    src = open(os.path.join(HERE, script), encoding="utf-8").read()
    for a, b in subs:
        assert a in src, (script, a)
        src = src.replace(a, b)
    exec(compile(src, script, "exec"), {"__name__": "__main__", "__file__": os.path.join(HERE, script)})


run("14_biomart.py", [(r"oz_tv_sic_work\genes14", r"oz_tv_sic_work\genes17")])
run("14_score_genes.py", [(r"oz_tv_sic_work\genes14", r"oz_tv_sic_work\genes17"),
                          ('os.path.join(T, "14_gene_intervals.csv")', 'os.path.join(T, "17_gene_intervals.csv")'),
                          ('os.path.join(T, "14_genes_all_intervals.csv")', 'os.path.join(W, "score_genes_all.csv")'),
                          ('os.path.join(T, "14_candidate_top3.csv")', 'os.path.join(W, "score_top3.csv")'),
                          ('os.path.join(T, "14_criterion_a_counts.csv")', 'os.path.join(W, "score_criterion_a_counts.csv")'),
                          ('os.path.join(T, "14_1B_dev_stress_genes.csv")', 'os.path.join(W, "score_1B_unused.csv")'),
                          ('IV[IV.locus.str.startswith("1B")].interval_id.iloc[0]', '"I06"')])
new = pd.read_csv(os.path.join(G, "score_genes_all.csv"))
old = pd.read_csv(os.path.join(OUT, "tables", "14_genes_all_intervals.csv"))
old = old[old.interval_id.isin(["I01", "I02", "I05", "I06", "I07", "I08"])]
A = pd.concat([old, new], ignore_index=True)
A.to_csv(os.path.join(G, "genes_all9.csv"), index=False, encoding="utf-8-sig")
print(A.groupby("interval_id").agg(n=("svevo_gene", "size"), n_a=("pts_a", lambda s: int((s > 0).sum()))))
