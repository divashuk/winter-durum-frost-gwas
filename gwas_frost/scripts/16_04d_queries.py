# -*- coding: utf-8 -*-
"""Раздел 16, шаг 4г. Запросы для поиска ортологов и копий по сборкам.
  * Гены-кандидаты (координаты Langdon):
      - локусы CBF Langdon (16_04_cbf_loci.py) в 5A 530–537 Мб, 5B 520–539 Мб и в повторе 5A 512–516 Мб,
        если идентичность лучшего белка ≥ 0.60 и покрытие ≥ 0.70 (фрагменты и псевдогены с идентичностью < 0.6 — только в таблице);
      - гены Svevo, перенесённые на Langdon в разделах 10/14 (genes10/genes_wide.csv, genes14/genes.csv; ид. ≥ 0.98, покр. ≥ 0.90):
        TRITD5Av1G171440, TRITD5Bv1G171760 (PRR95), TRITD7Bv1G107490, TRITD7Bv1G098240, TRITD1Bv1G205280,
        TRITD5Bv1G171580 и TRITD5Bv1G171610 (ближайшие гены к 5B:535723924 / 535618213).
    Запрос = ген + 2 кб вверх + 1 кб вниз по цепи гена (5′ → 3′ гена).
  * Якоря синтении: до 6 ближайших перенесённых генов Svevo (тело гена ≥ 1 кб) с каждой стороны группы в пределах 1 Мб (если их меньше 3 — 2 Мб)
    (тело гена Langdon); сторона L — проксимальнее, R — дистальнее по координате Langdon.
Выход: WSL s16/q/candidates.tsv, anchors.tsv, regions.txt (для samtools faidx) и scripts для извлечения."""
import pandas as pd, numpy as np, sys

sys.stdout.reconfigure(encoding="utf-8")
WK = "/work/scratch/"
S16 = WK + "s16/"
import os; os.makedirs(S16 + "q", exist_ok=True)

L = pd.read_csv(S16 + "langdon_cbf/langdon_cbf_loci.tsv", sep="\t")
inwin = ((L.chr == "5A") & L.start.between(530e6, 537e6)) | ((L.chr == "5B") & L.start.between(520e6, 539e6)) | \
        ((L.chr == "5A") & L.start.between(512e6, 516e6))
L = L[inwin].copy()
L["is_gene"] = (L.best_ident >= 0.6) & (L.best_cov >= 0.7)
L["cand_id"] = "CBF_" + L.chr + "_" + (L.start / 1e6).map(lambda x: f"{x:.3f}")
L["group"] = np.where((L.chr == "5A") & (L.start < 516e6), "5A_513dup", np.where(L.chr == "5A", "FrA2", "FrB2"))
L.to_csv(S16 + "q/langdon_cbf_candidates_all.tsv", sep="\t", index=False)

g10 = pd.read_csv(WK + "genes10/genes_wide.csv")[["chr", "langdon_start", "langdon_end", "strand", "svevo_gene", "identity", "coverage"]]
g14 = pd.read_csv(WK + "genes14/genes.csv")[["chr", "langdon_start", "langdon_end", "strand", "svevo_gene", "identity", "coverage"]]
G = pd.concat([g10, g14]).drop_duplicates("svevo_gene").sort_values(["chr", "langdon_start"]).reset_index(drop=True)

other = {"TRITD5Av1G171440": "5A_171440", "TRITD5Bv1G171760": "5B_535", "TRITD5Bv1G171580": "5B_535", "TRITD5Bv1G171610": "5B_535",
         "TRITD7Bv1G107490": "7B_348", "TRITD7Bv1G098240": "7B_282", "TRITD1Bv1G205280": "1B_653"}
# TRITD5Bv1G171430 (vWA, CS TraesCS5B02G319900) — в 9.8 кб от 5B:535618213; в разделе 14 не перенесён (идентичность кДНК 0.886 < 0.98);
# положение — minimap2 -x splice:hq кДНК TRITD5Bv1G171430.1 на Langdon 5B:535.3–536.0 Мб (s16/snp535)
G = pd.concat([G, pd.DataFrame([dict(chr="5B", langdon_start=535628037, langdon_end=535630040, strand="+", svevo_gene="TRITD5Bv1G171430",
                                     identity=0.886, coverage=1.0)])], ignore_index=True)
other["TRITD5Bv1G171430"] = "5B_535"
# вторая копия того же семейства: miniprot белка TRITD5Bv1G171430 на Langdon 5B:535.3–536.0 (s16/snp535) — 535719364–535721589 (+),
# идентичность 0.826, без сдвига рамки; в 2.3 кб от 5B:535723924. Ищется по синтенным окнам второго прохода (16_04l_local_pass2.sh)
G = pd.concat([G, pd.DataFrame([dict(chr="5B", langdon_start=535719364, langdon_end=535721589, strand="+", svevo_gene="TRITD5Bv1G171430_c2",
                                     identity=0.826, coverage=1.0)])], ignore_index=True)
other["TRITD5Bv1G171430_c2"] = "5B_535"
O = G[G.svevo_gene.isin(other)].copy()
O["group"] = O.svevo_gene.map(other)
cand = pd.concat([
    pd.DataFrame(dict(cand_id=L.cand_id, group=L.group, chr=L.chr, start=L.start, end=L.end, strand=L.strand, source="miniprot (Langdon)",
                      is_gene=L.is_gene, cs_paralog=L.cs_paralog, tm_name=L.tm_name, complete=L.complete)),
    pd.DataFrame(dict(cand_id=O.svevo_gene, group=O.group, chr=O.chr, start=O.langdon_start, end=O.langdon_end, strand=O.strand,
                      source="перенос Svevo v1 (разделы 10/14)", is_gene=True, cs_paralog="", tm_name="", complete=""))])
cand["qstart"] = np.where(cand.strand == "+", cand.start - 2000, cand.start - 1000)
cand["qend"] = np.where(cand.strand == "+", cand.end + 1000, cand.end + 2000)
cand.to_csv(S16 + "q/candidates.tsv", sep="\t", index=False, lineterminator="\n")

# якоря: по группам
anc = []
for grp, x in cand[cand.is_gene].groupby("group"):
    c, lo, hi = x.chr.iloc[0], x.start.min(), x.end.max()
    g = G[(G.chr == c) & (G.identity >= 0.98) & ~G.svevo_gene.isin(other) & (G.langdon_end - G.langdon_start >= 1000)]
    left = g[(g.langdon_end < lo) & (g.langdon_end > lo - 1e6)].tail(6)
    if len(left) < 3: left = g[(g.langdon_end < lo) & (g.langdon_end > lo - 2e6)].tail(6)       # мало генов — окно 2 Мб
    right = g[(g.langdon_start > hi) & (g.langdon_start < hi + 1e6)].head(6)
    if len(right) < 3: right = g[(g.langdon_start > hi) & (g.langdon_start < hi + 2e6)].head(6)
    for side, y in [("L", left), ("R", right)]:
        for _, r in y.iterrows():
            anc.append(dict(group=grp, anchor=r.svevo_gene, side=side, chr=c, start=int(r.langdon_start), end=int(r.langdon_end), strand=r.strand))
A = pd.DataFrame(anc).drop_duplicates(["group", "anchor"])
A.to_csv(S16 + "q/anchors.tsv", sep="\t", index=False, lineterminator="\n")
print(cand.groupby(["group", "is_gene"]).size()); print(A.groupby(["group", "side"]).size())
with open(S16 + "q/extract.sh", "w", newline="\n") as f:
    f.write("set -e; P=/work/scratch/primers11; Q=/work/scratch/s16/q; E=/opt/micromamba/envs/gwas_frost/bin/n: > $Q/queries.fa/n")
    for _, r in cand[cand.is_gene].iterrows():
        rc = "-i " if r.strand == "-" else ""
        f.write(f"$E/samtools faidx {rc}$P/langdon_{r.chr}.fa langdon_{r.chr}:{int(r.qstart)}-{int(r.qend)} | sed '1s/.*/>Q|{r.cand_id}|{r.group}/' >> $Q/queries.fa\n")
    for _, r in A.iterrows():
        f.write(f"$E/samtools faidx $P/langdon_{r.chr}.fa langdon_{r.chr}:{r.start}-{r.end} | sed '1s/.*/>A|{r.anchor}|{r.group}|{r.side}/' >> $Q/queries.fa\n")
    f.write("grep -c '>' $Q/queries.fa\n")
