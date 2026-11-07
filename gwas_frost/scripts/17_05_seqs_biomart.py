# -*- coding: utf-8 -*-
"""Раздел 17, шаг 4а (повтор 16.3 для генов части A). Ensembl Plants BioMart (функция mart из 16_03_biomart_cbf.py) для генов
tables/17_candidate_genes_other_loci.csv: biomart/svevo_interpro.tsv (InterPro с координатами на белке), svevo_cs_homologs.tsv,
svevo_positions.tsv; cs_interpro.tsv — для ортолога CS.
Последовательности: Ensembl REST /lookup и /sequence/id в момент работы отвечали HTTP 500 (/sequence/region работал), поэтому вместо
16_03b_rest_svevo_seqs.py — файлы Ensembl Plants 63 того же выпуска:
  Svevo — FTP pep.all и cds.all (WSL ref/, скачаны в разделе 17) и cdna.all (ref/, раздел 5); канонический транскрипт — метка
  tag=Ensembl_canonical в GFF3 (ref/), как is_canonical в REST;
  CS — Triticum_aestivum.IWGSC.pep.all.fa.gz на сервере (~/fgum_rasshirenie_2026/ref, источник белков CS в разделе 16); транскрипт .1, иначе первый.
Выход: WSL ~/oz_tv_sic_work/s17/{biomart,seq}/."""
import os, sys, re, gzip, subprocess
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
S17 = "/work/scratch/s17/"
REF = "/work/scratch/ref/"
os.makedirs(S17 + "biomart", exist_ok=True); os.makedirs(S17 + "seq", exist_ok=True)


def load(script, stop):
    """функции из скрипта раздела 16 (код до строки stop, без выполнения основной части)"""
    src = open(os.path.join(HERE, script), encoding="utf-8").read()
    ns = {"__name__": "lib"}
    exec(compile(src[:src.index(stop)], script, "exec"), ns)
    return ns


C = pd.read_csv(os.path.join(OUT, "tables/17_candidate_genes_other_loci.csv"))
genes = sorted(set(C.ген))
cs = sorted(set(C.ID_CS_RefSeq_v1_1.dropna()) - {""})
W = S17 + "biomart/"
if not os.path.exists(W + "cs_interpro.tsv"):
    bm = load("16_03_biomart_cbf.py", "def save(")
    save = lambda d, n: (d.to_csv(W + n + ".tsv", sep="\t", index=False), print(n, len(d)))[0]
    save(bm["mart"]("tturgidum_eg_gene", {"ensembl_gene_id": ",".join(genes)},
                    ["ensembl_gene_id", "taestivum_eg_homolog_ensembl_gene", "taestivum_eg_homolog_orthology_type", "taestivum_eg_homolog_perc_id",
                     "taestivum_eg_homolog_perc_id_r1", "taestivum_eg_homolog_chromosome", "taestivum_eg_homolog_chrom_start"]), "svevo_cs_homologs")
    save(bm["mart"]("tturgidum_eg_gene", {"ensembl_gene_id": ",".join(genes)}, ["ensembl_gene_id", "interpro", "interpro_description", "interpro_start",
                                                                              "interpro_end"]), "svevo_interpro")
    save(bm["mart"]("tturgidum_eg_gene", {"ensembl_gene_id": ",".join(genes)}, ["ensembl_gene_id", "chromosome_name", "start_position", "end_position",
                                                                              "strand", "description"]), "svevo_positions")
    save(bm["mart"]("taestivum_eg_gene", {"ensembl_gene_id": ",".join(cs)}, ["ensembl_gene_id", "interpro", "interpro_description", "interpro_start",
                                                                           "interpro_end"]), "cs_interpro")


def rdfa(path, keep):
    d, n = {}, None
    for line in gzip.open(path, "rt"):
        if line.startswith(">"):
            t = line[1:].split()[0]; n = t if t in keep else None
            if n: d[n] = []
        elif n: d[n].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}


can, rows = {}, []
for chrom in sorted({g[5:7] for g in genes}):
    for line in gzip.open(REF + f"Triticum_turgidum.Svevo.v1.63.chromosome.{chrom}.gff3.gz", "rt"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "mRNA" and "Ensembl_canonical" in f[8]:
            t = re.search(r"transcript_id=([^;]+)", f[8]).group(1); g = t.rsplit(".", 1)[0]
            if g in genes:
                can[g] = t; rows.append(dict(gene=g, transcript=t, chr=f[0], start=int(f[3]), end=int(f[4]), strand=1 if f[6] == "+" else -1))
keep = set(can.values())
seqs = {k: rdfa(REF + f"Triticum_turgidum.Svevo.v1.{k}.all.fa.gz", keep) for k in ["pep", "cds", "cdna"]}
for k, v in seqs.items():
    open(S17 + f"seq/svevo_{k}.fa", "w", newline="\n").write("".join(f">{g} {t}\n{v[t]}\n" for g, t in can.items() if t in v))
pd.DataFrame(rows).to_csv(S17 + "seq/svevo_canonical.tsv", sep="\t", index=False, lineterminator="\n")
print("Svevo: генов", len(can), {k: len(v) for k, v in seqs.items()})
pat = "|".join(cs)
cmd = (f"ssh compute-host 'zcat ~/fgum_rasshirenie_2026/ref/Triticum_aestivum.IWGSC.pep.all.fa.gz' | "
       f"/opt/micromamba/envs/gwas_frost/bin/seqkit grep -n -r -p 'gene:({pat}) '")
txt = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "--", "bash", "-c", cmd], capture_output=True, text=True).stdout
best = {}
for r in [r for r in txt.split(">") if r.strip()]:
    h, *sq = r.split("\n"); g = re.search(r"gene:(\S+)", h).group(1); t = h.split()[0]
    if g not in best or t.endswith(".1"): best[g] = (h, "".join(sq))
open(S17 + "seq/cs_pep.fa", "w", newline="\n").write("".join(f">{h}\n{s}\n" for h, s in best.values()))
print("CS белков", len(best), "из", len(cs), sorted(set(cs) - set(best)))
