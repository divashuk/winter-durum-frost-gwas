# -*- coding: utf-8 -*-
"""Этап 10.4. Сводка аннотации и балльная оценка генов широких интервалов.
Критерии (баллы; источник каждого утверждения — ID InterPro/GO, ID ортолога Compara или DOI, проверенный в Crossref):
 C1 GO холода (2): GO:0009409, GO:0009631, GO:0050826, GO:0070417 (BioMart Ensembl Plants, tturgidum_eg_gene)
 C2 CBF/DREB1 (3): InterPro «Dehydration-responsive element-binding protein 1»; прочие AP2/ERF — 1
 C3 ICE (3): ортолог Arabidopsis ICE1/SCRM2
 C4 CAMTA (3): ортолог CAMTA1–6 или InterPro CG-1 / CAMTA
 C5 COR/LEA/дегидрины (2): InterPro Dehydrin / Late embryogenesis abundant; ортологи COR15A/B, COR47, ERD10, ERD14, RD29A/B, KIN1/2, COR413, LTI
 C6 десатуразы жирных кислот (2): InterPro «Fatty acid desaturase» / ацил-ACP-десатураза (стерол- и каротин-десатуразы — нет); ортологи FAD2–FAD8, SSI2/FAB2, ADS
 C7 сахара/пролин (2): InterPro P5CS, раффинозосинтаза, трегалоза, сахарозосинтаза, SPS, GH32 (инвертазы/фруктаны); ортологи P5CS, GolS, RS, SUS, SPS, TPS, TPP, CWINV, VIN, BAM3
 C8 антиоксиданты (1): InterPro SOD, каталаза, аскорбатпероксидаза, глутатионредуктаза/-пероксидаза, DHAR, MDAR, пероксиредоксин
 C9 фитохромы/VRN/FT/PPD/CO-like (2; MADS без ортолога AP1/FUL/SOC1/SVP — 1): InterPro Phytochrome, CCT, B-box, PEBP, MADS; ортологи PHY, CO/COL/BBX, FT/TSF, AP1/FUL/CAL, PRR/TOC1, ELF3, GI, LHY/CCA1, SOC1, SVP, PIF
     (имена VRN1/VRN2 Arabidopsis не учитываются: это не гомологи VRN-1/VRN-2 злаков)
 C10 ортолог с опубликованной ролью в холодостойкости (2): список COLD_ORTHO с DOI; гены CBF на Fr-A2/Fr-B2 и VRN-1 — по DOI работ на пшенице
 C11 холодовая экспрессия ортолога CS в expVIP (1): SRP043554, log2FC >= 1 и p (Уэлч) < 0.05
 P1 в узком интервале (1); P2 ведущий SNP в гене или <= 100 т.п.н. (1); P3 SNP GBS в гене с −log10 p MLM/GEMMA >= Meff (1)
Кандидат — ген хотя бы с одним функциональным критерием C1–C11."""
import os, re, sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
W = r"/work/scratch/genes10"
T = r"/work/project/gwas_frost/tables"
rd = lambda n: pd.read_csv(os.path.join(W, n), sep="\t", dtype=str)
g = pd.read_csv(os.path.join(W, "genes_wide.csv"))
L = pd.read_csv(os.path.join(T, "10_loci.csv"))
ipr, go = rd("biomart_interpro.tsv"), rd("biomart_go.tsv")
hom = {s: rd(f"biomart_hom_{s}.tsv") for s in ["athaliana", "osativa", "taestivum"]}
desc = {s: rd(f"biomart_desc_{s}.tsv") for s in ["athaliana", "osativa"]}
csb = rd("cs_best_ortholog.tsv")
E = rd("expvip_SRP043554.tsv") if os.path.exists(os.path.join(W, "expvip_SRP043554.tsv")) else pd.DataFrame(columns=["cs_gene"])

DOI = {"Knox2010": "10.1007/s00122-010-1288-7", "Pearce2013": "10.1007/s00122-013-2165-y", "Sieber2016": "10.1007/s00122-016-2685-3",
       "Wurschum2017": "10.1111/tpj.13424", "Dhillon2010": "10.1104/pp.110.159079", "Li2015_expVIP_SRP043554": "10.1105/tpc.114.134338",
       "Close1996": "10.1034/j.1399-3054.1996.970422.x"}
COLD_ORTHO = {
    "AT4G25490": ("CBF1/DREB1B", "10.1126/science.280.5360.104"), "AT3G26744": ("ICE1", "10.1101/gad.1077503"),
    "AT5G09410": ("CAMTA1", "10.1105/tpc.108.063958"), "AT5G64220": ("CAMTA2", "10.1105/tpc.108.063958"), "AT2G22300": ("CAMTA3", "10.1105/tpc.108.063958"),
    "AT1G56070": ("LOS1", "10.1073/pnas.112040099"), "AT3G53110": ("LOS4", "10.1073/pnas.172399299"), "AT3G55990": ("ESK1", "10.1111/j.1365-313X.2006.02994.x"),
    "AT5G05580": ("FAD8", "10.1104/pp.106.4.1615"), "AT3G12120": ("FAD2", "10.1073/pnas.90.13.6208"),
    "AT2G47180": ("GolS1", "10.1046/j.0960-7412.2001.01227.x"), "AT1G56600": ("GolS2", "10.1046/j.0960-7412.2001.01227.x"), "AT1G09350": ("GolS3", "10.1046/j.0960-7412.2001.01227.x"),
    "AT5G24470": ("PRR5", "10.1093/pcp/pcp004"), "AT5G02810": ("PRR7", "10.1093/pcp/pcp004"), "AT2G46790": ("PRR9", "10.1093/pcp/pcp004"),
    "AT2G43010": ("PIF4", "10.1073/pnas.1211295109"), "AT2G18790": ("PHYB", "10.1038/ng.2007.3")}
KNOWN = {"VRN-A1": ("5A", 598268462, 598273527), "VRN-B1": ("5B", 605129460, 605142609),
         "Fr-A2 (CBF)": ("5A", 533000000, 535200000), "Fr-B2 (CBF)": ("5B", 523000000, 524200000)}

# --- сводные строки аннотации на ген
ipr = ipr.dropna(subset=["interpro"])
IPR = ipr.groupby("ensembl_gene_id").apply(lambda d: "; ".join(sorted(set(d.interpro + " " + d.interpro_description.fillna(""))))).rename("interpro")
GO = go.dropna(subset=["go_id"]).groupby("ensembl_gene_id").apply(lambda d: "; ".join(sorted(set(d.go_id + " " + d.name_1006.fillna(""))))).rename("go")


def homstr(sp, n=3):
    h = hom[sp].dropna(subset=[f"{sp}_eg_homolog_ensembl_gene"]).copy()
    h["pid"] = h[f"{sp}_eg_homolog_perc_id"].astype(float)
    if sp in desc:
        h = h.merge(desc[sp].rename(columns={"ensembl_gene_id": f"{sp}_eg_homolog_ensembl_gene", "external_gene_name": "nm", "description": "ds"}), how="left")
    else:
        h["nm"] = h[f"{sp}_eg_homolog_associated_gene_name"]
    h["nm"] = h["nm"].fillna(h[f"{sp}_eg_homolog_associated_gene_name"]).fillna("")
    h = h.sort_values("pid", ascending=False)
    s = h.groupby("ensembl_gene_id").apply(lambda d: "; ".join(f"{a} {b} ({c.replace('ortholog_', '')}, {p:.0f}%)" for a, b, c, p in
                                         list(zip(d[f"{sp}_eg_homolog_ensembl_gene"], d.nm, d[f"{sp}_eg_homolog_orthology_type"], d.pid))[:n]))
    return h, s


hA, A = homstr("athaliana"); hO, O = homstr("osativa")
csb = csb.rename(columns={"ensembl_gene_id": "svevo_gene", "taestivum_eg_homolog_ensembl_gene": "cs_gene", "taestivum_eg_homolog_orthology_type": "cs_type",
                          "pid": "cs_pid", "taestivum_eg_homolog_chromosome": "cs_chr", "taestivum_eg_homolog_chrom_start": "cs_start"})
csb["cs_type"] = csb.cs_type.str.replace("ortholog_", "")
G = g.merge(csb[["svevo_gene", "cs_gene", "cs_type", "cs_pid", "cs_chr", "cs_start"]], how="left", on="svevo_gene")
G["interpro"] = G.svevo_gene.map(IPR).fillna(""); G["go"] = G.svevo_gene.map(GO).fillna("")
G["at_orthologs"] = G.svevo_gene.map(A).fillna(""); G["os_orthologs"] = G.svevo_gene.map(O).fillna("")
if len(E):
    E = E.astype({"tpm_cold": float, "tpm_control": float, "log2fc_cold": float, "p_welch": float})
    G = G.merge(E[["cs_gene", "tpm_cold", "tpm_control", "log2fc_cold", "p_welch"]], how="left", on="cs_gene")
else:
    G[["tpm_cold", "tpm_control", "log2fc_cold", "p_welch"]] = np.nan
at_names = hA.groupby("ensembl_gene_id").apply(lambda d: list(zip(d.athaliana_eg_homolog_ensembl_gene, d.nm))).to_dict()
thr = {r.locus_id: (4.54 if r.ld_set == "S1" else 4.87) for r in L.itertuples()}


def has(s, pat):
    m = re.findall(r"(IPR\d+) ([^;]*(?:%s)[^;]*)" % pat, s, flags=re.I)
    return ", ".join(sorted({a for a, _ in m}))


def score(r):
    ev, pts = [], {}
    ats = at_names.get(r.svevo_gene, [])
    atn = {n.upper(): a for a, n in ats if isinstance(n, str)}
    def at_match(pat):
        return ", ".join(f"{a} {n}" for a, n in ats if isinstance(n, str) and re.match(pat, n, flags=re.I))
    # C1
    cg = [x.split(" ")[0] for x in r.go.split("; ") if x.split(" ")[0] in ("GO:0009409", "GO:0009631", "GO:0050826", "GO:0070417")]
    if cg: pts["C1"] = 2; ev.append("C1 GO холода: " + ", ".join(cg))
    # C2
    d = has(r.interpro, "Dehydration-responsive element-binding protein 1")
    if d: pts["C2"] = 3; ev.append(f"C2 CBF/DREB1: {d}")
    else:
        a = has(r.interpro, "AP2/ERF domain")
        if a: pts["C2"] = 1; ev.append(f"C2 AP2/ERF (не DREB1): {a}")
    # C3
    m = at_match(r"^(ICE\d|SCRM\d?)$")
    if m: pts["C3"] = 3; ev.append("C3 ICE: ортолог " + m)
    # C4
    m = at_match(r"^CAMTA\d$"); d = has(r.interpro, "CG-1|CAMTA|Calmodulin-binding transcription")
    if m or d: pts["C4"] = 3; ev.append("C4 CAMTA: " + ", ".join(x for x in [d, m] if x))
    # C5
    d = has(r.interpro, r"Dehydrin|Late embryogenesis abundant|\bLEA\b|LEA domain"); m = at_match(r"^(COR\d+\w*|ERD1[04]|RD29[AB]|KIN\d|COR413\w*|LTI\d+|XERO\d|HIRD11)$")
    if d or m: pts["C5"] = 2; ev.append("C5 COR/LEA/дегидрин: " + ", ".join(x for x in [d, m] if x) + (f" (дегидрины: doi:{DOI['Close1996']})" if "ehydrin" in r.interpro else ""))
    # C6
    d = has(r.interpro, r"Fatty acid desaturase|Acyl-\[acyl-carrier-protein\] desaturase|Acyl-ACP desaturase"); m = at_match(r"^(FAD\d|SSI2|FAB2|ADS\d)$")   # стерол- и каротин-десатуразы не учитываются
    if d or m: pts["C6"] = 2; ev.append("C6 десатураза: " + ", ".join(x for x in [d, m] if x))
    # C7
    d = has(r.interpro, "pyrroline-5-carboxylate|Raffinose synthase|Trehalose(?! [Tt]ransporter)|Sucrose synthase|Sucrose-phosphate|Glycoside hydrolase family 32|Glycosyl hydrolase family 32|Fructan|Galactinol")
    m = at_match(r"^(P5CS\d|P5CR|GolS\d|RS\d|SUS\d|SPS\w*|TPS\d+|TPP[A-J]?|CWINV\d|VIN\d|BAM3)$")
    if d or m: pts["C7"] = 2; ev.append("C7 сахара/пролин: " + ", ".join(x for x in [d, m] if x))
    # C8
    d = has(r.interpro, "Superoxide dismutase|Catalase|scorbate peroxidase|Glutathione reductase|Glutathione peroxidase|ehydroascorbate reductase|Peroxiredoxin")
    m = at_match(r"^(CSD\d|FSD\d|MSD\d|CAT\d|APX\d|GR\d?|GPX\d|DHAR\d|MDAR\d)$")
    if d or m: pts["C8"] = 1; ev.append("C8 антиоксидант: " + ", ".join(x for x in [d, m] if x))
    # C9
    m = at_match(r"^(PHY[A-E]|CO|COL\d+|BBX\d+|FT|TSF|AP1|AGL8|FUL|CAL|TOC1|ELF\d|GI|LHY|CCA1|SOC1|SVP|PIF\d)$")
    # PRR — только псевдорегуляторы ответа по AGI (имена PRR1/PRR2 в Arabidopsis носят и пинорезинолредуктазы AT1G32100/AT4G13660)
    prr = ", ".join(f"{a} {n}" for a, n in ats if a in ("AT5G60100", "AT5G24470", "AT5G02810", "AT2G46790", "AT5G61380"))
    m = ", ".join(x for x in [m, prr] if x)
    d = has(r.interpro, "Phytochrome|CCT domain|B-box|Phosphatidylethanolamine-binding")
    dm = has(r.interpro, "MADS")
    if m or d: pts["C9"] = 2; ev.append("C9 свет/цветение: " + ", ".join(x for x in [d, dm, m] if x))
    elif dm: pts["C9"] = 1; ev.append("C9 MADS-box: " + dm)
    # C10
    c10 = [f"{a} {COLD_ORTHO[a][0]} doi:{COLD_ORTHO[a][1]}" for a, _ in ats if a in COLD_ORTHO]
    for k, (ch, s, e) in KNOWN.items():
        near = r.chr == ch and r.langdon_end >= s - 5e4 and r.langdon_start <= e + 5e4
        if near and k.startswith("Fr-A2") and "C2" in pts and pts["C2"] == 3:
            c10.append(f"ген CBF в кластере {k}: doi:{DOI['Knox2010']}; doi:{DOI['Sieber2016']}; doi:{DOI['Wurschum2017']}")
        if near and k.startswith("Fr-B2") and "C2" in pts and pts["C2"] == 3:
            c10.append(f"ген CBF в кластере {k}: doi:{DOI['Pearce2013']}")
        if near and k.startswith("VRN") and at_match(r"^(AP1|AGL8|FUL|CAL)$"):
            c10.append(f"{k} (ортолог AP1/FUL на позиции гена): doi:{DOI['Dhillon2010']}")
    if c10: pts["C10"] = 2; ev.append("C10 ортолог/ген с ролью в холодостойкости: " + " | ".join(dict.fromkeys(c10)))
    # C11
    if pd.notna(r.log2fc_cold):
        if r.log2fc_cold >= 1 and r.p_welch < 0.05:
            pts["C11"] = 1; ev.append(f"C11 expVIP SRP043554: {r.cs_gene} log2FC={r.log2fc_cold}, p={r.p_welch:.3g} (doi:{DOI['Li2015_expVIP_SRP043554']})")
        expn = f"{r.log2fc_cold} (p={r.p_welch:.2g})"
    else:
        expn = "нет данных expVIP" if pd.notna(r.cs_gene) else "нет ортолога CS"
    func = sum(pts.values())
    # близость
    pr = {}
    if r.in_narrow: pr["P1"] = 1
    if r.dist_to_lead <= 1e5: pr["P2"] = 1
    if pd.notna(r.gbs_best_mlp_in_gene) and r.gbs_best_mlp_in_gene >= thr[r.locus_id]: pr["P3"] = 1
    return pd.Series({"func_score": func, "prox_score": sum(pr.values()), "total_score": func + sum(pr.values()),
                      "criteria": ",".join(list(pts) + list(pr)), "evidence": " || ".join(ev), "expvip_log2fc_cold": expn})


G = pd.concat([G, G.apply(score, axis=1)], axis=1)
G["dist_to_lead_kb"] = (G.dist_to_lead / 1e3).round(1)
cols = ["locus_id", "locus_name", "category", "chr", "langdon_start", "langdon_end", "strand", "svevo_gene", "svevo_chr", "svevo_start", "biotype",
        "identity", "coverage", "in_narrow", "dist_to_lead_kb", "gbs_snps_in_gene", "gbs_best_mlp_in_gene", "gbs_best_at", "cs_gene", "cs_type", "cs_pid",
        "cs_chr", "cs_start", "at_orthologs", "os_orthologs", "interpro", "go", "tpm_cold", "tpm_control", "expvip_log2fc_cold",
        "func_score", "prox_score", "total_score", "criteria", "evidence"]
G = G[cols].sort_values(["locus_id", "langdon_start"])
G.to_csv(os.path.join(T, "10_genes_wide.csv"), index=False, encoding="utf-8-sig")
G[G.in_narrow].to_csv(os.path.join(T, "10_genes_narrow.csv"), index=False, encoding="utf-8-sig")
C = G[G.func_score > 0].sort_values(["locus_id", "total_score", "dist_to_lead_kb"], ascending=[True, False, True])
C.to_csv(os.path.join(T, "10_candidate_genes.csv"), index=False, encoding="utf-8-sig")
S = G.groupby(["locus_id", "locus_name"]).agg(n_wide=("svevo_gene", "size"), n_narrow=("in_narrow", "sum"), n_cand_wide=("func_score", lambda s: (s > 0).sum()),
                                              n_cand_narrow=("func_score", lambda s: 0), expvip=("tpm_cold", lambda s: s.notna().sum())).reset_index()
S["n_cand_narrow"] = S.locus_id.map(G[G.in_narrow & (G.func_score > 0)].groupby("locus_id").size()).fillna(0).astype(int)
print(S.to_string())
print("кандидатов всего:", len(C), "; в узких:", int(C.in_narrow.sum()))
