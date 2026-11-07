# -*- coding: utf-8 -*-
"""Раздел 17, шаг 7и. Проверка ведущих SNP 1B, 7B и 2B на паралоги -> tables/17_lead_snp_paralog_check.csv.
  * число попаданий шаблона Langdon ± 150 п.н. (BLAST шаблонов 17_14b, evalue 1e-10) с идентичностью ≥ 90 % в Langdon (14 хромосом) и
    CS RefSeq v1.0 (весь геном): всего и с выравниванием ≥ 150 п.н. (≥ половины шаблона); целевое попадание Langdon — в том числе;
  * доля гетерозигот по SNP у 326 образцов потока 2 (исходный VCF коллекции, GT) и средняя по всем 35 406 SNP VCF у тех же образцов
    (17_11b_parent_gt.sh: s17/kasp/het_by_snp.tsv, het_genome.tsv). Для сравнения приведены и остальные SNP типа (б)."""
import os, sys, importlib
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
K = lib.S17 + "kasp/"
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qseq", "sseq"]
H = {g: pd.read_csv(K + f"tmpl_{g}.tsv", sep="\t", header=None, names=COLS) for g in ["langdon", "cs"]}
het = pd.read_csv(K + "het_by_snp.tsv", sep="\t", header=None, names=["chr", "pos", "n_called", "n_het", "het"])
gm = pd.read_csv(K + "het_genome.tsv", sep="\t", header=None).iloc[0, 1]
SB = pd.read_csv(K + "snps_b.tsv", sep="\t")
rows = []
for _, s in SB[~SB.role.str.startswith("заместитель")].iterrows():
    mk = f"KASP_{s.chr}_{int(s.pos)}"; r = dict(SNP=s.snp, локус=s.locus, проверка_паралогов="да (задание)" if s.chr in ("1B", "7B", "2B") else "для сравнения")
    for g, d in H.items():
        x = d[(d.qseqid == mk) & (d.pident >= 90)]
        r[f"попаданий_{g}_ид90"] = len(x); r[f"попаданий_{g}_ид90_150пн"] = int((x.length >= 150).sum())
        r[f"хромосомы_{g}"] = ", ".join(sorted(set(x[x.length >= 150].sseqid.astype(str).str.replace("langdon_", "", regex=False))))
    h = het[(het.chr == s.chr) & (het.pos == int(s.pos))]
    r.update(образцов_с_вызовом=int(h.n_called.iloc[0]) if len(h) else None, гетерозигот=int(h.n_het.iloc[0]) if len(h) else None,
             доля_гетерозигот=round(float(h.het.iloc[0]), 4) if len(h) else None, средняя_доля_гетерозигот_по_VCF=round(float(gm), 4))
    r["отношение_к_средней"] = round(r["доля_гетерозигот"] / gm, 2) if r["доля_гетерозигот"] is not None else None
    rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(os.path.join(lib.OUT, "tables/17_lead_snp_paralog_check.csv"), index=False, encoding="utf-8-sig")
print(T.to_string())
