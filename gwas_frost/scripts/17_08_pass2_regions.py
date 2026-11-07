# -*- coding: utf-8 -*-
"""Раздел 17, шаг 4д. Синтенные окна, копии и участки второго прохода — scripts/16_04g_pass2_regions.py (правила 16_04f_copies.py) без
изменения метода, рабочая папка s17. PAF: s17/scan/<сборка>.paf (сервер, 17_07a_server_scan.sh), Svevo — s17/svevo/svevo.paf.
Запуск: python scripts/17_08_pass2_regions.py  (все сборки, для которых есть PAF)."""
import os, sys, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
GEN = ["Langdon", "Svevo", "Kronos", "ASM5281831v1", "Chili", "Mahmoudi", "Zavitan_WEW2.1", "Dicoccum_lpTriDico3", "CS_v1.0", "Jagger", "Julius",
       "SY_Mattis", "Norin61", "ArinaLrFor", "Stanley", "Landmark"]
args = []
for g in GEN:
    p = lib.S17 + ("svevo/svevo.paf" if g == "Svevo" else f"scan/{g}.paf")
    if os.path.exists(p): args.append(f"{g}={p}")
lib.run16("16_04g_pass2_regions.py", [("S16 = cp.S16", f'S16 = "{lib.S17}"')], argv=args)
