# -*- coding: utf-8 -*-
"""Раздел 17, шаг 7д. Подбор праймеров KASP — scripts/16_07d_kasp_design.py без изменения метода (Tm 58–62 °C, ампликон 50–120 п.н.,
маска, специфичность к гомеологам и паралогам по BLAST шаблонов, ΔG > −9 ккал/моль), рабочая папка s17/kasp.
Повтор для 5B:535618213 (задание 17): тот же код с Tm специфичной части 57–63 °C и ампликоном 50–150 п.н. (designs_535.tsv);
в итоговый designs.tsv этот маркер берётся из повтора."""
import os, sys, importlib
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
K = lib.S17 + "kasp/"
S = ('S16 = "/work/scratch/s16/"', f'S16 = "{lib.S17}"')
MK = "KASP_5B_535618213"
lib.run16("16_07d_kasp_design.py", [S, ('T = T[T.status == "мишень выбрана"].copy()', f'T = T[(T.status == "мишень выбрана") & (T.marker != "{MK}")].copy()')])
os.replace(K + "designs.tsv", K + "designs_main.tsv"); os.replace(K + "primers_cand.fa", K + "primers_cand_main.fa")
lib.run16("16_07d_kasp_design.py", [S, ('T = T[T.status == "мишень выбрана"].copy()', f'T = T[(T.status == "мишень выбрана") & (T.marker == "{MK}")].copy()'),
                                    ("if not (58 <= tx <= 62 and 58 <= ty <= 62)", "if not (57 <= tx <= 63 and 57 <= ty <= 63)"),
                                    ("if not (50 <= amp <= 120): continue", "if not (50 <= amp <= 150): continue"),
                                    ("if not 58 <= tc <= 62", "if not 57 <= tc <= 63"),
                                    ('K + "designs.tsv"', 'K + "designs_535.tsv"'), ('K + "primers_cand.fa"', 'K + "primers_cand_535.fa"')])
D = pd.concat([pd.read_csv(K + "designs_main.tsv", sep="\t"), pd.read_csv(K + "designs_535.tsv", sep="\t")], ignore_index=True)
D.to_csv(K + "designs.tsv", sep="\t", index=False, lineterminator="\n")
open(K + "primers_cand.fa", "w", newline="\n").write(open(K + "primers_cand_main.fa").read() + open(K + "primers_cand_535.fa").read())
print(D[D.marker == MK][["marker", "status"]].head())
