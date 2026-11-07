# -*- coding: utf-8 -*-
"""Этап 10.6 (в). Q1B.6 (L08): гены с доменами LRR, NB-ARC, киназными, Rx/CC N-концевыми, TIR — состав доменов по InterPro.
Вывод о принадлежности к NLR не делается: перечисляется только, какие домены есть у гена."""
import os, re, sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
T = r"/work/project/gwas_frost/tables"
G = pd.read_csv(os.path.join(T, "10_genes_wide.csv"))
g = G[G.locus_id == "L08"].copy()
DOM = {"LRR": r"Leucine-rich repeat", "NB-ARC": r"NB-ARC", "P-loop NTPase": r"P-loop containing nucleoside triphosphate hydrolase",
       "Rx N-концевой (CC)": r"Rx.*N-terminal|Disease resistance protein, Rx|Rx, N-terminal", "TIR": r"Toll/interleukin-1 receptor|\bTIR\b",
       "киназный домен": r"Protein kinase domain|kinase-like domain", "WRKY": r"WRKY", "Disease resistance (описание InterPro)": r"Disease resistance"}
for k, p in DOM.items():
    g[k] = g.interpro.fillna("").apply(lambda s: ", ".join(sorted({a for a, d in re.findall(r"(IPR\d+) ([^;]*)", s) if re.search(p, d, re.I)})))
sel = g[(g["LRR"] != "") | (g["NB-ARC"] != "") | (g["Disease resistance (описание InterPro)"] != "")].copy()
sel["состав доменов"] = sel.apply(lambda r: " + ".join(k for k in DOM if r[k]), axis=1)
cols = ["svevo_gene", "chr", "langdon_start", "langdon_end", "in_narrow", "dist_to_lead_kb", "gbs_snps_in_gene", "gbs_best_mlp_in_gene", "cs_gene",
        "at_orthologs", "os_orthologs", "состав доменов"] + list(DOM) + ["expvip_log2fc_cold"]
sel[cols].sort_values("langdon_start").to_csv(os.path.join(T, "10_Q1B6_LRR_genes.csv"), index=False, encoding="utf-8-sig")
print("генов в Q1B.6:", len(g), "; с LRR/NB-ARC/«Disease resistance»:", len(sel), "; из них в узком интервале:", int(sel.in_narrow.sum()))
print(sel[["svevo_gene", "langdon_start", "in_narrow", "dist_to_lead_kb", "состав доменов", "gbs_best_mlp_in_gene"]].sort_values("langdon_start").to_string())
