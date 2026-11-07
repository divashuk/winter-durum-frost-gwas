# -*- coding: utf-8 -*-
"""Раздел 17, шаг 4б (повтор 16.4, правила 16_04d_queries.py). Запросы для поиска копий генов части A по 16 сборкам.
  * гены — tables/17_candidate_genes_other_loci.csv, кроме TRITD1Bv1G205280 (разобран в разделе 16, результаты берутся оттуда);
    запрос = ген + 2 кб вверх + 1 кб вниз по цепи гена; группа синтении — своя у каждого гена;
  * якоря синтении — до 6 ближайших перенесённых генов Svevo (тело ≥ 1 кб, ид. ≥ 0.98) с каждой стороны в пределах 1 Мб (если < 3 — 2 Мб);
    источник переноса — разделы 10, 14, 17 (genes10/genes_wide.csv, genes14/genes.csv, genes17/genes.csv);
  * белки для miniprot (второй проход): белок канонического транскрипта Svevo (SV|<ген>) и ортолога CS (CS|<ген CS>) — s17/seq;
  * Svevo v1: в разделе 16 Svevo был только 5A, 5B и участки ±20 кб; здесь для синтении — участки Svevo от крайнего якоря до крайнего
    ± 1.1 Мб (Ensembl REST /sequence/region, ≤ 10 Мб на запрос; целые хромосомы не скачиваются); для 5B — из ref11/svevo_5B.fa.
Выход: WSL s17/q/candidates.tsv, candidate_names.tsv, anchors.tsv, extract.sh, prot_long.fa, prot2.fa (пустой: CBF нет), svevo_windows.tsv."""
import os, re, gzip, sys
import pandas as pd, numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
WK = "/work/scratch/"
S17 = WK + "s17/"; os.makedirs(S17 + "q", exist_ok=True)
C = pd.read_csv(os.path.join(OUT, "tables/17_candidate_genes_other_loci.csv"))
C = C[C.ген != "TRITD1Bv1G205280"].copy()
m = C.позиция_Langdon.str.extract(r"^(\w+):(\d+)-(\d+) \(([+-])\)$")
C["chr"], C["start"], C["end"], C["strand"] = m[0], m[1].astype(int), m[2].astype(int), m[3]
C["group"] = C.chr + "_" + (C.start / 1e6).astype(int).astype(str)
cand = pd.DataFrame(dict(cand_id=C.ген, group=C.group, chr=C.chr, start=C.start, end=C.end, strand=C.strand,
                         source="перенос Svevo v1 (разделы 14/17)", is_gene=True, cs_paralog="", tm_name="", complete=""))
cand["qstart"] = np.where(cand.strand == "+", cand.start - 2000, cand.start - 1000)
cand["qend"] = np.where(cand.strand == "+", cand.end + 1000, cand.end + 2000)
cand.to_csv(S17 + "q/candidates.tsv", sep="\t", index=False, lineterminator="\n")
names = pd.DataFrame(dict(локус=C.локус, ген=C.ген, cand_id=C.ген, ID_Svevo=C.ген, ID_CS_RefSeq_v1_1=C.ID_CS_RefSeq_v1_1.fillna(""), CS_идентичность="",
                          позиция_Langdon=C.позиция_Langdon, функция=C.функция, основание=C.обоснование,
                          приоритет=np.where(C.роль == "основной", "1 (основной)", "2 (" + C.роль + ")")))
names.to_csv(S17 + "q/candidate_names.tsv", sep="\t", index=False, lineterminator="\n")
cols = ["chr", "langdon_start", "langdon_end", "strand", "svevo_gene", "identity", "coverage"]
G = pd.concat([pd.read_csv(WK + "genes10/genes_wide.csv")[cols], pd.read_csv(WK + "genes14/genes.csv")[cols],
               pd.read_csv(WK + "genes17/genes.csv")[cols]]).drop_duplicates("svevo_gene").sort_values(["chr", "langdon_start"]).reset_index(drop=True)
anc = []
for _, x in cand.iterrows():
    c, lo, hi = x.chr, x.start, x.end
    g = G[(G.chr == c) & (G.identity >= 0.98) & ~G.svevo_gene.isin(set(cand.cand_id)) & (G.langdon_end - G.langdon_start >= 1000)]
    left = g[(g.langdon_end < lo) & (g.langdon_end > lo - 1e6)].tail(6)
    if len(left) < 3: left = g[(g.langdon_end < lo) & (g.langdon_end > lo - 2e6)].tail(6)
    right = g[(g.langdon_start > hi) & (g.langdon_start < hi + 1e6)].head(6)
    if len(right) < 3: right = g[(g.langdon_start > hi) & (g.langdon_start < hi + 2e6)].head(6)
    for side, y in [("L", left), ("R", right)]:
        for _, r in y.iterrows():
            anc.append(dict(group=x.group, anchor=r.svevo_gene, side=side, chr=c, start=int(r.langdon_start), end=int(r.langdon_end), strand=r.strand))
A = pd.DataFrame(anc).drop_duplicates(["group", "anchor"])
A.to_csv(S17 + "q/anchors.tsv", sep="\t", index=False, lineterminator="\n")
print(A.groupby(["group", "side"]).size().unstack())
with open(S17 + "q/extract.sh", "w", newline="\n") as f:
    f.write("set -e; P=/work/scratch/primers11; Q=/work/scratch/s17/q; E=/opt/micromamba/envs/gwas_frost/bin/n: > $Q/queries.fa/n")
    for _, r in cand.iterrows():
        rc = "-i " if r.strand == "-" else ""
        f.write(f"$E/samtools faidx {rc}$P/langdon_{r.chr}.fa langdon_{r.chr}:{int(r.qstart)}-{int(r.qend)} | sed '1s/.*/>Q|{r.cand_id}|{r.group}/' >> $Q/queries.fa\n")
    for _, r in A.iterrows():
        f.write(f"$E/samtools faidx $P/langdon_{r.chr}.fa langdon_{r.chr}:{r.start}-{r.end} | sed '1s/.*/>A|{r.anchor}|{r.group}|{r.side}/' >> $Q/queries.fa\n")
    f.write("grep '>' $Q/queries.fa | grep '^>Q' | sed 's/>//' > $Q/q_ids.txt; $E/seqkit grep -f $Q/q_ids.txt $Q/queries.fa > $Q/queries_Q.fa; grep -c '>' $Q/queries.fa $Q/queries_Q.fa\n")


def fa(path):
    d, n = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"): n = line[1:].strip(); d[n] = []
        else: d[n].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}


sv = {k.split()[0]: v for k, v in fa(S17 + "seq/svevo_pep.fa").items()}
cs = {re.search(r"gene:(\S+)", k).group(1): v for k, v in fa(S17 + "seq/cs_pep.fa").items()}
with open(S17 + "q/prot_long.fa", "w", newline="\n") as f:
    for _, r in C.iterrows():
        f.write(f">SV|{r.ген}\n{sv[r.ген].rstrip('*')}\n")
        if isinstance(r.ID_CS_RefSeq_v1_1, str) and r.ID_CS_RefSeq_v1_1 in cs:
            f.write(f">CS|{r.ID_CS_RefSeq_v1_1}\n{cs[r.ID_CS_RefSeq_v1_1].rstrip('*')}\n")
open(S17 + "q/prot2.fa", "w").close()
# участки Svevo: положения генов Svevo по GFF3 Ensembl Plants 63 (ref/)
pos = {}
for c in sorted(set(cand.chr)):
    with gzip.open(WK + f"ref/Triticum_turgidum.Svevo.v1.63.chromosome.{c}.gff3.gz", "rt") as fh:
        for line in fh:
            f_ = line.split("\t")
            if len(f_) > 8 and f_[2] == "gene":
                gid = re.search(r"ID=gene:([^;]+)", f_[8]).group(1); pos[gid] = (f_[0], int(f_[3]), int(f_[4]))
win = []
for grp, a in A.groupby("group"):
    ids = list(a.anchor) + list(cand[cand.group == grp].cand_id)
    p = [pos[i] for i in ids if i in pos]
    ch = pd.Series([x[0] for x in p]).mode()[0]
    s = [x for x in p if x[0] == ch]
    med = np.median([x[1] for x in s]); s = [x for x in s if abs(x[1] - med) <= 5e6]
    lo, hi = max(1, int(min(x[1] for x in s) - 1.1e6)), int(max(x[2] for x in s) + 1.1e6)
    win.append(dict(group=grp, chr=ch, start=lo, end=hi, n_svevo_genes=len(s), length_Mb=round((hi - lo) / 1e6, 2)))
W = pd.DataFrame(win)
W.to_csv(S17 + "q/svevo_windows.tsv", sep="\t", index=False, lineterminator="\n")
print(W.to_string())
