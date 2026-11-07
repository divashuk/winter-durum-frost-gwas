# -*- coding: utf-8 -*-
"""Раздел 17, шаг 7д′. Второй тур подбора для мишеней, у которых в 17_14c праймеры не подобраны (SNP типа б и гены без годного маркера):
тот же код 16_07d с условиями повтора 5B:535618213 (Tm специфичной части 57–63 °C, ампликон 50–150 п.н.). В итоговый designs.tsv для этих
мишеней берутся варианты второго тура (designs.tsv первого тура сохраняется как designs_round1.tsv). Праймеры второго тура —
primers_cand_relaxed.fa (BLAST — 17_14b с Q=primers_cand_relaxed.fa, результаты дописываются к primer_*.tsv)."""
import os, sys, importlib
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
K = lib.S17 + "kasp/"
D = pd.read_csv(K + "designs.tsv", sep="\t"); T = pd.read_csv(K + "targets.tsv", sep="\t")
ok_m = set(D[D.status == "подобран"].marker)
P = pd.read_csv(os.path.join(lib.OUT, "tables/17_kasp_primers.csv"))
chosen_genes = set(P[P.маркер.notna() & P.праймер_FAM.notna()].ген.dropna())
T = T[T.status == "мишень выбрана"]
need = T[(~T.marker.isin(ok_m)) & ((T.type.str.startswith("б")) | (~T.gene.isin(chosen_genes) & T.gene.notna()))].marker.tolist()
need = [m for m in need if m != "KASP_5B_535618213"]
print("второй тур:", len(need), need)
open(K + "relaxed_markers.txt", "w").write("\n".join(need) + "\n")
lib.run16("16_07d_kasp_design.py", [
    ('S16 = "/work/scratch/s16/"', f'S16 = "{lib.S17}"'),
    ('T = T[T.status == "мишень выбрана"].copy()', f'T = T[(T.status == "мишень выбрана") & T.marker.isin({need!r})].copy()'),
    ("if not (58 <= tx <= 62 and 58 <= ty <= 62)", "if not (57 <= tx <= 63 and 57 <= ty <= 63)"),
    ("if not (50 <= amp <= 120): continue", "if not (50 <= amp <= 150): continue"),
    ("if not 58 <= tc <= 62", "if not 57 <= tc <= 63"),
    ('K + "designs.tsv"', 'K + "designs_relaxed.tsv"'), ('K + "primers_cand.fa"', 'K + "primers_cand_relaxed.fa"')])
R = pd.read_csv(K + "designs_relaxed.tsv", sep="\t")
if not os.path.exists(K + "designs_round1.tsv"): D.to_csv(K + "designs_round1.tsv", sep="\t", index=False, lineterminator="\n")
D = pd.read_csv(K + "designs_round1.tsv", sep="\t")
M = pd.concat([D[~D.marker.isin(need)], R.assign(note="второй тур: Tm 57–63 °C, ампликон ≤ 150 п.н.")], ignore_index=True)
M.to_csv(K + "designs.tsv", sep="\t", index=False, lineterminator="\n")
print(R.groupby("marker").status.first().to_string())
