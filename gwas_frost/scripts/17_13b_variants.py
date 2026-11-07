# -*- coding: utf-8 -*-
"""Раздел 17, шаги 5б–6. Варианты, классы функциональности и группы гаплотипов генов части A — scripts/16_06_variants.py без изменения
метода (замены путей — 17_lib.VAR_SUBS): tables/17_variants_functional.csv, tables/17_haplotype_groups.csv.
TRITD1Bv1G205280 (1B) разобран в разделе 16: его строки из tables/16_variants_functional.csv и 16_haplotype_groups.csv дописываются
с пометкой источника."""
import os, sys, importlib
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
lib.run16("16_06_variants.py", lib.VAR_SUBS)
for t in ["variants_functional", "haplotype_groups"]:
    a = pd.read_csv(os.path.join(lib.OUT, f"tables/17_{t}.csv")).assign(источник="раздел 17")
    b = pd.read_csv(os.path.join(lib.OUT, f"tables/16_{t}.csv"))
    b = b[b.cand_id == "TRITD1Bv1G205280"].assign(источник=f"раздел 16 (tables/16_{t}.csv)")
    pd.concat([a, b], ignore_index=True).to_csv(os.path.join(lib.OUT, f"tables/17_{t}.csv"), index=False, encoding="utf-8-sig")
