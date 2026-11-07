# -*- coding: utf-8 -*-
"""Раздел 17, шаг 6а. Ортологи генов части A в ячмене, рисе и Brachypodium (Ensembl Compara, REST /homology) —
scripts/16_06a_compara.py без изменения метода; рабочая папка s17 (s17/compara/orthologs.tsv, orthologs_pep.fa).
Все гены части A — гены Svevo (ID TRITD…), таблица локусов CBF Langdon не нужна."""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
lib.run16("16_06a_compara.py", [('S16 = "/work/scratch/s16/"', f'S16 = "{lib.S17}"'),
                                (r'L = pd.read_csv(S16 + "langdon_cbf/langdon_cbf_loci.tsv", sep="\t")', "L = None")])
