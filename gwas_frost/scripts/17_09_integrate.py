# -*- coding: utf-8 -*-
"""Раздел 17, шаг 4з. Сводка копий генов части A по 16 сборкам и последовательности — scripts/16_04k_integrate.py без изменения метода
(копии по ДНК 16_04f, белок/CDS по miniprot, «ORF полная», гомеолог D у мягкой, неразмещённые контиги). Замены: рабочая папка s17;
genes16 -> genes17; tables/16_gene_copies_by_genome.csv -> tables/17_gene_copies_by_genome.csv; правило «Svevo 7B/1B — участок ±20 кб,
синтения не проверялась» не применяется: в 17 у Svevo участки с якорями (17_07b_svevo.sh), синтения проверяется как у прочих сборок.
TRITD1Bv1G205280 (1B) разобран в разделе 16: его строки из tables/16_gene_copies_by_genome.csv дописываются в таблицу 17 с пометкой."""
import os, sys, importlib
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
lib.run16("16_04k_integrate.py", [
    ("S16 = cl.WSL", f'S16 = "{lib.S17}"'),
    ('os.path.join(OUT, "genes16")', 'os.path.join(OUT, "genes17")'),
    ('os.path.join(OUT, "genes16", f"{base}_{kind}.fa")', 'os.path.join(OUT, "genes17", f"{base}_{kind}.fa")'),
    ('os.path.join(OUT, "tables/16_cbf_loci_by_genome.csv")', f'"{lib.S17}copies/cbf_loci_unused.csv"'),
    ('os.path.join(OUT, "tables/16_gene_copies_by_genome.csv")', 'os.path.join(OUT, "tables/17_gene_copies_by_genome.csv")'),
    ('if g == "Svevo" and r.tchr in ("7B", "1B"):', 'if False:'),
    ('d = cl.parse_gff(S16 + f"p2/{g}_win.gff", plen, chr_map=lambda s: s)', 'd = pd.DataFrame()   # белков CBF в разделе 17 нет')])
T = pd.read_csv(os.path.join(lib.OUT, "tables/17_gene_copies_by_genome.csv"))
T["источник"] = "раздел 17"
T16 = pd.read_csv(os.path.join(lib.OUT, "tables/16_gene_copies_by_genome.csv"))
T16 = T16[T16.cand_id == "TRITD1Bv1G205280"].assign(источник="раздел 16 (tables/16_gene_copies_by_genome.csv)")
pd.concat([T, T16], ignore_index=True).to_csv(os.path.join(lib.OUT, "tables/17_gene_copies_by_genome.csv"), index=False, encoding="utf-8-sig")
