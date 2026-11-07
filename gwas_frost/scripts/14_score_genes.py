# -*- coding: utf-8 -*-
"""Этап 14.4–14.5. Аннотация и приоритизация генов интервалов 14.
Описание функции гена — только из Ensembl Plants (BioMart, выпуск 63): InterPro (ID + описание), GO (ID), ортологи Compara в Arabidopsis
(AGI, имя, описание TAIR из athaliana_eg_gene) и рисе; ортолог IWGSC RefSeq v1.1 (гены TraesCS…02G в Ensembl Plants) — лучший по % идентичности,
с приоритетом той же хромосомы. Собственных догадок нет.
Баллы:
 (а) 3 — связь с холодом / морозостойкостью / яровизацией по описанию: CBF/DREB1 (IPR045277; ортолог DDF1/DDF2/CBF1–4), COR/LEA (InterPro «Late embryogenesis
     abundant»; ортологи COR15A/B, COR47, ERD10, ERD14, RD29A/B, KIN1/2, COR413), дегидрины (IPR000167), ICE (ICE1/SCRM2), CAMTA (IPR005559; CAMTA1–6),
     VRN (ортолог AP1/FUL/AGL8/CAL), PPD (PRR3/5/7/9, TOC1 по AGI), фитохромы (InterPro «Phytochrome»; PHYA–E), PIF (PIF1–8, PIL), FT (PEBP, IPR008914),
     CO/COL/BBX (CCT IPR010402 или B-box IPR000315 с ортологом CO/COL/BBX), GO холода (GO:0009409, GO:0009631, GO:0050826, GO:0070417);
 (б) близость: ведущий SNP в гене или <= 100 т.п.н. — 2; <= 1 Мб — 1;
 (в) LD: r² ближайшего к гену SNP с ведущим >= 0.8 — 2; >= 0.5 — 1.
Итог = (а) + (б) + (в). Топ-3 на интервал (при равенстве — ближе к ведущему SNP)."""
import os, re, sys
import pandas as pd
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
W = r"/work/scratch/genes14"
T = r"/work/project/gwas_frost/tables"
rd = lambda n: pd.read_csv(os.path.join(W, n), sep="\t", dtype=str)
g = pd.read_csv(os.path.join(W, "genes.csv"))
IV = pd.read_csv(os.path.join(T, "14_gene_intervals.csv"))
ipr, go = rd("biomart_interpro.tsv").dropna(subset=["interpro"]), rd("biomart_go.tsv").dropna(subset=["go_id"])
hA, hO, hT = rd("biomart_hom_athaliana.tsv"), rd("biomart_hom_osativa.tsv"), rd("biomart_hom_taestivum.tsv")
dA, dO, dT = rd("biomart_desc_athaliana.tsv"), rd("biomart_desc_osativa.tsv"), rd("biomart_desc_taestivum.tsv")
IPR = ipr.groupby("ensembl_gene_id").apply(lambda d: "; ".join(sorted(set(d.interpro + " " + d.interpro_description.fillna(""))))).to_dict()
GO = go.groupby("ensembl_gene_id").go_id.apply(lambda s: ",".join(sorted(set(s)))).to_dict()
hA = hA.dropna(subset=["athaliana_eg_homolog_ensembl_gene"]).merge(dA.rename(columns={"ensembl_gene_id": "athaliana_eg_homolog_ensembl_gene", "external_gene_name": "at_name", "description": "at_desc"}), how="left")
hA["pid"] = hA.athaliana_eg_homolog_perc_id.astype(float); hA = hA.sort_values("pid", ascending=False)
hA["at_name"] = hA.at_name.fillna(hA.athaliana_eg_homolog_associated_gene_name).fillna("")
hA["at_desc"] = hA.at_desc.fillna("").str.replace(r"\s*\[Source:.*\]", "", regex=True)
AT = hA.groupby("ensembl_gene_id").apply(lambda d: list(zip(d.athaliana_eg_homolog_ensembl_gene, d.at_name, d.at_desc, d.athaliana_eg_homolog_orthology_type.str.replace("ortholog_", ""), d.pid))[:3]).to_dict()
hO = hO.dropna(subset=["osativa_eg_homolog_ensembl_gene"]).merge(dO.rename(columns={"ensembl_gene_id": "osativa_eg_homolog_ensembl_gene", "external_gene_name": "os_name", "description": "os_desc"}), how="left")
hO["pid"] = hO.osativa_eg_homolog_perc_id.astype(float); hO = hO.sort_values("pid", ascending=False)
OS = hO.groupby("ensembl_gene_id").apply(lambda d: "; ".join(f"{a} {('' if pd.isna(b) else b)}".strip() for a, b in list(zip(d.osativa_eg_homolog_ensembl_gene, d.os_name))[:2])).to_dict()
hT = hT.dropna(subset=["taestivum_eg_homolog_ensembl_gene"]).copy(); hT["pid"] = hT.taestivum_eg_homolog_perc_id.astype(float)
chrmap = g.drop_duplicates("svevo_gene").set_index("svevo_gene").chr.to_dict()
hT["same"] = [c == chrmap.get(s, "") for s, c in zip(hT.ensembl_gene_id, hT.taestivum_eg_homolog_chromosome)]
CS = hT.sort_values(["same", "pid"], ascending=False).drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
csdesc = dT.set_index("ensembl_gene_id").description.fillna("").str.replace(r"\s*\[Source:.*\]", "", regex=True).to_dict()
PRR = {"AT5G60100": "PRR3", "AT5G24470": "PRR5", "AT5G02810": "PRR7", "AT2G46790": "PRR9", "AT5G61380": "TOC1"}
COLD_GO = {"GO:0009409", "GO:0009631", "GO:0050826", "GO:0070417"}


def cat_a(gid):
    s = IPR.get(gid, ""); ats = AT.get(gid, []); names = [(a, n) for a, n, *_ in ats]; hits = []
    def ip(pat):
        return sorted({a for a, d in re.findall(r"(IPR\d+) ([^;]*)", s) if re.search(pat, d, re.I)})
    def at(pat):
        return [f"{a} {n}" for a, n in names if isinstance(n, str) and re.match(pat, n, re.I)]
    x = ip(r"Dehydration-responsive element-binding protein 1") + at(r"^(DDF[12]|CBF\d|DREB1\w?)$")
    if x: hits.append("CBF/DREB1: " + ", ".join(x))
    x = ip(r"Late embryogenesis abundant") + at(r"^(COR\d+\w*|ERD1[04]|RD29[AB]|KIN\d|COR413\w*)$")
    if x: hits.append("COR/LEA: " + ", ".join(x))
    x = ip(r"Dehydrin")
    if x: hits.append("дегидрин: " + ", ".join(x))
    x = at(r"^(ICE\d|SCRM\d?)$")
    if x: hits.append("ICE: " + ", ".join(x))
    x = ip(r"CG-1|CAMTA") + at(r"^CAMTA\d$")
    if x: hits.append("CAMTA: " + ", ".join(x))
    mads = ip(r"MADS"); x = at(r"^(AP1|AGL8|FUL|CAL)$")
    if x and mads: hits.append("VRN-1-подобный (AP1/FUL): " + ", ".join(mads + x))
    x = [f"{a} {PRR[a]}" for a, _ in names if a in PRR]
    if x: hits.append("PPD/PRR: " + ", ".join(x))
    x = ip(r"Phytochrome") + at(r"^PHY[A-E]$")
    if x: hits.append("фитохром: " + ", ".join(x))
    x = at(r"^(PIF\d|PIL\d)$")
    if x: hits.append("PIF: " + ", ".join(x))
    x = ip(r"Phosphatidylethanolamine-binding protein$|PEBP")
    if x: hits.append("FT-подобный (PEBP): " + ", ".join(x))
    x = at(r"^(CO|COL\d+|BBX\d+)$")
    if x and ip(r"CCT domain|B-box"): hits.append("CO/COL/BBX: " + ", ".join(ip(r"CCT domain|B-box") + x))
    cg = sorted(set(GO.get(gid, "").split(",")) & COLD_GO)
    if cg: hits.append("GO холода: " + ", ".join(cg))
    return " | ".join(hits)


def desc(gid):
    parts = []
    if IPR.get(gid): parts.append("InterPro: " + IPR[gid])
    for a, n, d, t, p in AT.get(gid, [])[:2]:
        parts.append(f"Arabidopsis {a} {n} ({t}, {p:.0f}%): {d}".strip())
    if OS.get(gid): parts.append("рис: " + OS[gid])
    return " || ".join(parts) if parts else "нет InterPro и ортологов в Ensembl Plants"


g["cs_refseq_v11"] = g.svevo_gene.map(lambda x: CS.taestivum_eg_homolog_ensembl_gene.get(x, ""))
g["cs_type"] = g.svevo_gene.map(lambda x: str(CS.taestivum_eg_homolog_orthology_type.get(x, "")).replace("ortholog_", ""))
g["cs_pid"] = g.svevo_gene.map(lambda x: CS.pid.get(x, np.nan))
g["cs_description"] = g.cs_refseq_v11.map(lambda x: csdesc.get(x, ""))
g["function_description"] = g.svevo_gene.map(desc)
g["source"] = "Ensembl Plants 63 BioMart: InterPro, GO, Compara (Arabidopsis, рис, T. aestivum IWGSC RefSeq v1.1)"
g["go_ids"] = g.svevo_gene.map(lambda x: GO.get(x, ""))
g["criterion_a"] = g.svevo_gene.map(cat_a)
g["pts_a"] = np.where(g.criterion_a != "", 3, 0)
g["pts_b"] = np.where(g.dist_to_lead <= 1e5, 2, np.where(g.dist_to_lead <= 1e6, 1, 0))
g["pts_c"] = np.where(g.r2_nearest_with_lead >= 0.8, 2, np.where(g.r2_nearest_with_lead >= 0.5, 1, 0))
g["score"] = g.pts_a + g.pts_b + g.pts_c
g["dist_to_lead_kb"] = (g.dist_to_lead / 1e3).round(1)
g = g.merge(IV[["interval_id", "lead_set", "status", "start_Mb", "end_Mb"]], on="interval_id")
cols = ["interval_id", "locus", "lead_set", "status", "chr", "start_Mb", "end_Mb", "lead", "svevo_gene", "langdon_start", "langdon_end", "strand", "svevo_pos", "cs_refseq_v11", "cs_type", "cs_pid",
        "cs_description", "function_description", "go_ids", "dist_to_lead_kb", "nearest_snp", "r2_nearest_with_lead", "criterion_a", "pts_a", "pts_b", "pts_c", "score", "source"]
g = g[cols].sort_values(["interval_id", "langdon_start"])
g.to_csv(os.path.join(T, "14_genes_all_intervals.csv"), index=False, encoding="utf-8-sig")
top = g.sort_values(["interval_id", "score", "dist_to_lead_kb"], ascending=[True, False, True]).groupby("interval_id").head(3)
top.to_csv(os.path.join(T, "14_candidate_top3.csv"), index=False, encoding="utf-8-sig")
S = g.groupby(["interval_id", "locus", "lead_set"]).agg(n_genes=("svevo_gene", "size"), n_criterion_a=("pts_a", lambda s: int((s > 0).sum()))).reset_index()
S.to_csv(os.path.join(T, "14_criterion_a_counts.csv"), index=False, encoding="utf-8-sig")
print(S.to_string())
for iid, d in top.groupby("interval_id"):
    print("==", iid, d.locus.iloc[0], d.lead_set.iloc[0])
    for _, r in d.iterrows():
        print(f"  {r.svevo_gene} {r.langdon_start/1e6:.3f} Мб d={r.dist_to_lead_kb} кб r2={r.r2_nearest_with_lead} score={r.score} ({r.pts_a}/{r.pts_b}/{r.pts_c}) CS={r.cs_refseq_v11} | {r.criterion_a[:120]} | {r.function_description[:150]}")
# 1B: гены развития и стресса по аннотации (широкий список ключевых слов)
b = g[g.interval_id == IV[IV.locus.str.startswith("1B")].interval_id.iloc[0]].copy()
kw = r"Heat shock|HSF|WRKY|NAC domain|Myb|bZIP|AP2/ERF|NB-ARC|Leucine-rich repeat|Disease resistance|Dehydrin|Late embryogenesis|CCT|B-box|MADS|Phytochrome|Pseudo-response|Response regulator|Zinc finger, C2H2|Calmodulin|Protein kinase"
b["dev_stress_terms"] = b.function_description.str.findall(kw, flags=re.I).map(lambda l: ", ".join(sorted(set(x for x in l))))
b = b[(b.dev_stress_terms != "") | (b.criterion_a != "")]
b.to_csv(os.path.join(T, "14_1B_dev_stress_genes.csv"), index=False, encoding="utf-8-sig")
print("1B: генов с терминами развития/стресса:", len(b))
