# -*- coding: utf-8 -*-
"""Раздел 16, шаг 3. Список генов-кандидатов -> tables/16_candidate_genes.csv.
Источники: локусы CBF Langdon (miniprot, s16/langdon_cbf/langdon_cbf_loci.tsv), перенос генов Svevo (разделы 10/14),
BioMart Ensembl Plants (s16/biomart: InterPro, гомологи Svevo -> CS), ампликоны праймеров Pearce et al. 2013 на Langdon
(tables/12_pearce_primers_langdon_genome_amplicons.csv), белки CBF T. monococcum (Miller et al. 2006).
Имя CBF: (1) по ампликону праймера Pearce, лежащему в локусе (± 1 кб) — для 5B; (2) иначе CBF-<геном><n>, если идентичность
с белком TmCBF<n> ≥ 0.90; (3) иначе «CBF (CS <ID>)». Имена CS из Ensembl (external_gene_name) приводятся отдельно.
На 5B имена по ампликонам Pearce имеют приоритет; если имя по Tm уже занято таким локусом — «<имя>-подобный».
Приоритет: 1 — CBF Fr-A2 и Fr-B2, PRR95, TRITD5Av1G171440, гены у SNP 5B:535618213 / 535723924 (±50 кб); 2 — 7B; 3 — 1B;
«—» — фрагменты CBF (покрытие белка < 0.7 или идентичность < 0.6), в поиске копий не используются."""
import sys
import pandas as pd, numpy as np

sys.stdout.reconfigure(encoding="utf-8")
S16 = "/work/scratch/s16/"
C = pd.read_csv(S16 + "q/langdon_cbf_candidates_all.tsv", sep="\t")
Q = pd.read_csv(S16 + "q/candidates.tsv", sep="\t")
amp = pd.read_csv("tables/12_pearce_primers_langdon_genome_amplicons.csv")
amp = amp[amp.target == True]
cs = pd.read_csv(S16 + "biomart/cs_IPR045277_5ABD.tsv", sep="\t")
csname = dict(zip(cs.ensembl_gene_id, cs.external_gene_name.fillna("")))
ipr_sv = pd.read_csv(S16 + "biomart/svevo_interpro.tsv", sep="\t").drop_duplicates(["ensembl_gene_id", "interpro"])
ipr_cs = pd.read_csv(S16 + "biomart/cs_interpro.tsv", sep="\t").drop_duplicates(["ensembl_gene_id", "interpro"])
hom = pd.read_csv(S16 + "biomart/svevo_cs_homologs.tsv", sep="\t")


def ipr_txt(d, g):
    x = d[d.ensembl_gene_id == g]
    return "; ".join(f"{a} {b}" for a, b in zip(x.interpro, x.interpro_description)) if len(x) else ""


def cs_same_chr(g, chrom):
    x = hom[(hom.ensembl_gene_id == g) & (hom.taestivum_eg_homolog_chromosome == chrom)]
    return x.sort_values("taestivum_eg_homolog_perc_id", ascending=False).taestivum_eg_homolog_ensembl_gene.iloc[0] if len(x) else ""


rows = []
pearce_names = set()
for _, r in C.iterrows():
    a = amp[(amp.chr == r.chr) & (amp.start <= r.end + 1000) & (amp.end >= r.start - 1000)]
    pearce_names |= {g for g in set(a.gene.str.replace(r"-pair \d", "", regex=True)) if g.startswith(f"CBF-{r.chr[-1]}")}
for _, r in C.iterrows():
    a = amp[(amp.chr == r.chr) & (amp.start <= r.end + 1000) & (amp.end >= r.start - 1000)]
    pearce = ", ".join(sorted(g for g in set(a.gene.str.replace(r"-pair \d", "", regex=True)) if g.startswith(f"CBF-{r.chr[-1]}")))
    gl = r.chr[-1]
    tm_nm = f"CBF-{gl}{str(r.tm_name).replace('CBF', '')}"
    if pearce and r.chr == "5B":
        name, how = pearce, "ампликон праймеров Pearce et al. 2013 в локусе Langdon"
    elif isinstance(r.tm_name, str) and r.tm_ident >= 0.90 and r.chr == "5B" and tm_nm in pearce_names:
        name, how = tm_nm + "-подобный", f"белок T. monococcum {r.tm_name} (ид. {r.tm_ident:.3f}); имя {tm_nm} занято локусом с ампликоном Pearce"
    elif isinstance(r.tm_name, str) and r.tm_ident >= 0.90:
        name, how = f"CBF-{gl}{r.tm_name.replace('CBF', '')}", f"белок T. monococcum {r.tm_name}, идентичность {r.tm_ident:.3f}"
    else:
        name, how = f"CBF (CS {r.cs_paralog})", "ближайший ген CS той же хромосомы"
    sv = r.svevo_best if isinstance(r.svevo_best, str) and r.svevo_ident >= 0.97 and f"TRITD{r.chr}" in r.svevo_best else ""   # только та же хромосома
    grp = {"FrA2": "5A Fr-A2", "FrB2": "5B Fr-B2", "5A_513dup": "5A, повтор участка Fr-A2 (513.2 Мб)"}[r.group]
    basis = (f"miniprot белков CBF на Langdon; лучший белок {r.best_target} (ид. {r.best_ident:.3f}, покрытие {r.best_cov:.2f}); "
             f"{'полная копия' if r.complete == 'полная' else 'неполная: ' + str(r.why)}; имя: {how}"
             + (f"; имя в Ensembl CS: {csname.get(r.cs_paralog)}" if csname.get(r.cs_paralog) else "")
             + (f"; праймеры Pearce: {pearce}" if pearce and r.chr != "5B" else ""))
    rows.append(dict(локус=grp, ген=name, cand_id=r.cand_id, ID_Svevo=sv or "нет в Svevo v1 (≥ 0.97)", ID_CS_RefSeq_v1_1=r.cs_paralog,
                     CS_идентичность=r.cs_ident, позиция_Langdon=f"{r.chr}:{int(r.start)}-{int(r.end)} ({r.strand})",
                     функция=ipr_txt(ipr_cs, r.cs_paralog) or "DREB1/CBF (IPR045277)", основание=basis,
                     приоритет="1" if r.is_gene and r.group != "5A_513dup" else ("2 (повтор вне кластера)" if r.is_gene else "— (фрагмент)")))
other = {"TRITD5Av1G171440": ("5A Fr-A2", "ведущий SNP RU 5A:514681016 внутри гена (раздел 14, I03, топ-1)", "1"),
         "TRITD5Bv1G171760": ("5B Fr-B2", "PRR95; CCT + receiver; ортолог AT5G24470 PRR5 / AT2G46790 PRR9; 434 т.п.н. от 5B:535723924, r² = 1.0 (раздел 14)", "1"),
         "TRITD5Bv1G171430": ("5B Fr-B2", "копия 1 семейства в Langdon (5B:535.628 Мб): splice:hq кДНК Svevo, ид. 0.886 (в разделе 14 не перенесён); "
                              "miniprot белка Svevo — ид. 0.836, сдвиг рамки; в 9.8 кб от 5B:535618213", "1"),
         "TRITD5Bv1G171430_c2": ("5B Fr-B2", "копия 2 семейства TRITD5Bv1G171430 в Langdon (5B:535.719 Мб): miniprot белка Svevo, ид. 0.826, без сдвига рамки; "
                                 "в 2.3 кб от 5B:535723924 (SNP за 3′-концом гена)", "1"),
         "TRITD5Bv1G171580": ("5B Fr-B2", "ближайший ген к 5B:535723924 (91 кб; в ±50 кб генов нет)", "1 (ближайший, вне ±50 кб)"),
         "TRITD5Bv1G171610": ("5B Fr-B2", "второй ближайший ген к 5B:535723924 (174 кб)", "3 (контекст)"),
         "TRITD7Bv1G107490": ("7B 118–463 Мб", "топ-1 интервала I05: ведущий SNP IT 7B:348683836 в 208 т.п.н. (раздел 14)", "2"),
         "TRITD7Bv1G098240": ("7B 118–463 Мб", "единственный ген по критерию (а) в узком окне I08 (LEA_2 IPR004864), 66.7 Мб от SNP (раздел 14)", "2"),
         "TRITD1Bv1G205280": ("1B 650–661 Мб", "ведущий SNP RU 1B:653803203 внутри гена (NB-ARC + LRR), без маркера (раздел 14)", "3")}
for g, (loc, basis, pr) in other.items():
    q = Q[Q.cand_id == g].iloc[0]
    g0 = g.replace("_c2", "")
    rows.append(dict(локус=loc, ген=g if g == g0 else f"{g0}-подобный (копия 2)", cand_id=g, ID_Svevo=g0 if g == g0 else f"{g0} (паралог; отдельного гена в аннотации Svevo v1 нет)",
                     ID_CS_RefSeq_v1_1=cs_same_chr(g0, q.chr), CS_идентичность=None,
                     позиция_Langdon=f"{q.chr}:{int(q.start)}-{int(q.end)} ({q.strand})", функция=ipr_txt(ipr_sv, g0), основание=basis, приоритет=pr))
D = pd.DataFrame(rows)
D["_k"] = D.локус.map({"5A Fr-A2": 0, "5A, повтор участка Fr-A2 (513.2 Мб)": 1, "5B Fr-B2": 2, "7B 118–463 Мб": 3, "1B 650–661 Мб": 4})
D["_p"] = D.позиция_Langdon.str.extract(r":(\d+)-")[0].astype(int)
D = D.sort_values(["_k", "_p"])
dup = D.duplicated(["локус", "ген"], keep=False)   # одинаковое имя у двух локусов одной группы -> добавить положение
D.loc[dup, "ген"] = D.loc[dup, "ген"] + " @" + (D.loc[dup, "_p"] / 1e6).map(lambda x: f"{x:.3f}") + " Мб"
D = D.drop(columns=["_k", "_p"])
D.to_csv("tables/16_candidate_genes.csv", index=False, encoding="utf-8-sig")
D.to_csv(S16 + "q/candidate_names.tsv", sep="\t", index=False, lineterminator="\n")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 40)
print(D[["локус", "ген", "ID_Svevo", "ID_CS_RefSeq_v1_1", "позиция_Langdon", "приоритет"]].to_string())
