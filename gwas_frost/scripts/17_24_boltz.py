# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 4 (WSL, окружение s17p, GPU). Комплекс белок–ДНК Boltz-2 (Passaro et al. 2025) для генов CBF, у выбранной пары
которых есть замена или индел в AP2/ERF (InterPro IPR001471) или в сигнатурах PKK…ETRHP / DSAWR (s17/prot/diffs.tsv, 17_20a), — для каждого
аллеля пары: белок (полный, цепь A) + двухцепочечная ДНК CRT/DRE 5′-ATTTACCGACATAAG-3′ (цепь B) и комплементарная 5′-CTTATGTCGGTAAAT-3′ (цепь C).
MSA белка — файл .a3m ColabFold того же аллеля (сервер MMseqs2; 17_22), при его отсутствии — --use_msa_server; --no_kernels (ядра cuequivariance 0.13 несовместимы с torch 2.7.1); остальные параметры Boltz по умолчанию (recycling 3, sampling 200, 1 образец диффузии), вывод PDB.
Метрики (confidence_*.json): ipTM, pTM, pLDDT комплекса; контакты белок–ДНК: пары остаток–нуклеотид с тяжёлыми атомами в ≤ 4.0 Å;
изменение контактов между аллелями — по позициям, сопоставленным выравниванием белков (как в 17_20a).
Выход: s17/prot/boltz/<cand_id>_<аллель>/…, s17/prot/boltz_contacts.tsv, строки «Boltz-2» в tables/17_protein_structure_comparison.csv."""
import os, re, json, glob, subprocess, sys
import numpy as np, pandas as pd

W = "/work/scratch/s17/prot/"; B = W + "boltz/"; os.makedirs(B, exist_ok=True)
OUT = "/work/project/gwas_frost/"
EXE = "/opt/micromamba/envs/s17p/bin/boltz"
DNA1, DNA2 = "ATTTACCGACATAAG", "CTTATGTCGGTAAAT"
SP = pd.read_csv(W + "selected_pairs.tsv", sep="\t"); DF = pd.read_csv(W + "diffs.tsv", sep="\t"); DM = pd.read_csv(W + "domains.tsv", sep="\t")
seq, n = {}, None
for line in open(W + "alleles.fa"):
    if line.startswith(">"): n = line[1:].strip(); seq[n] = ""
    else: seq[n] += line.strip()
sel = []
for _, p in SP[SP.cand_id.str.startswith("CBF_")].iterrows():
    d = DF[(DF.cand_id == p.cand_id) & (DF.a == p.a) & (DF.b == p.b)]
    hit = d[d.домен.fillna("").str.contains("IPR001471|сигнатура")]
    if len(hit): sel.append((p, "; ".join(f"{r.aa_a}{int(r.pos_a) if not pd.isna(r.pos_a) else ''}{r.aa_b if r.тип == 'замена' else ' (' + r.тип + ')'}" for _, r in hit.iterrows())))
print("генов CBF для Boltz-2:", len(sel), flush=True)
for p, _ in sel:
    for al in (p.a, p.b):
        name = f"{p.cand_id}_{al}".replace(".", "_")
        if glob.glob(B + f"out/boltz_results_{name}/predictions/{name}/confidence_*.json"): continue
        y = B + f"{name}.yaml"
        a3m = W + f"cf/out/{p.cand_id}_{al}.a3m"            # MSA ColabFold того же аллеля (сервер MMseqs2)
        msa = f"      msa: {a3m}\n" if os.path.exists(a3m) else ""
        open(y, "w").write(f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {seq[f'{p.cand_id}|{al}']}\n{msa}"
                           f"  - dna:\n      id: B\n      sequence: {DNA1}\n  - dna:\n      id: C\n      sequence: {DNA2}\n")
        r = subprocess.run([EXE, "predict", y, "--out_dir", B + "out", *([] if msa else ["--use_msa_server"]), "--output_format", "pdb", "--override", "--no_kernels"],
                           capture_output=True, text=True)
        print(name, "код", r.returncode, r.stderr[-300:] if r.returncode else "", flush=True)


def atoms(pdb):
    out = []
    for l in open(pdb):
        if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H":
            out.append((l[21], int(l[22:26]), l[17:20].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


rows, crow = [], []
from Bio import Align
from Bio.Align import substitution_matrices
AL = Align.PairwiseAligner(); AL.substitution_matrix = substitution_matrices.load("BLOSUM62"); AL.open_gap_score, AL.extend_gap_score, AL.mode = -10, -0.5, "global"
for p, what in sel:
    res = {}
    for al in (p.a, p.b):
        name = f"{p.cand_id}_{al}".replace(".", "_")
        d = B + f"out/boltz_results_{name}/predictions/{name}/"
        js = sorted(glob.glob(d + "confidence_*model_0.json")); pdb = sorted(glob.glob(d + "*model_0.pdb"))
        if not js or not pdb: res[al] = None; continue
        c = json.load(open(js[0])); X = atoms(pdb[0])
        P_ = [x for x in X if x[0] == "A"]; N_ = [x for x in X if x[0] in "BC"]
        pc = np.array([x[3] for x in P_]); nc = np.array([x[3] for x in N_])
        dist = np.linalg.norm(pc[:, None] - nc[None], axis=2)
        cont = set()
        for i, j in zip(*np.where(dist <= 4.0)):
            cont.add((P_[i][1], P_[i][2], f"{N_[j][0]}:{N_[j][2]}{N_[j][1]}"))
        resid = sorted({(r, n_) for r, n_, _ in cont})
        res[al] = dict(conf=c, cont=cont, resid=resid)
        for r, n_ in resid:
            crow.append(dict(cand_id=p.cand_id, allele=al, pos=r, aa=n_, nucleotides=", ".join(sorted(x[2] for x in cont if x[0] == r))))
    if not all(res.get(al) for al in (p.a, p.b)):
        rows.append(dict(ген=p.ген, cand_id=p.cand_id, аллели=f"{p.a} vs {p.b}", метод="Boltz-2 (белок–ДНК CRT/DRE)", строка="комплекс", участок=what,
                         примечание="модель не получена")); continue
    aln = AL.align(seq[f"{p.cand_id}|{p.a}"], seq[f"{p.cand_id}|{p.b}"])[0]; m = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1): m[s1 + k + 1] = s2 + k + 1
    ra = {r for r, _ in res[p.a]["resid"]}; rb = {r for r, _ in res[p.b]["resid"]}
    only_a = sorted(r for r in ra if m.get(r) not in rb); only_b = sorted(r for r in rb if r not in {m.get(x) for x in ra})
    sa, sb = seq[f"{p.cand_id}|{p.a}"], seq[f"{p.cand_id}|{p.b}"]
    for al, oth in ((p.a, only_a), (p.b, only_b)):
        c = res[al]["conf"]; s = sa if al == p.a else sb
        rows.append(dict(ген=p.ген, cand_id=p.cand_id, аллели=f"{p.a} vs {p.b}", метод="Boltz-2 (белок–ДНК CRT/DRE)", строка="комплекс", участок=f"аллель {al}; замены в AP2/сигнатурах: {what}",
                         ipTM=c.get("iptm"), pTM=c.get("ptm"), pLDDT_комплекса=c.get("complex_plddt"),
                         контактов_остаток_нуклеотид=len(res[al]["cont"]), остатков_в_контакте=len(res[al]["resid"]),
                         состав_контактов=", ".join(f"{n_}{r}" for r, n_ in res[al]["resid"]),
                         контакты_только_у_этого_аллеля=", ".join(f"{s[r - 1]}{r}" for r in oth)))
    print(p.cand_id, "готово", flush=True)
pd.DataFrame(crow).to_csv(W + "boltz_contacts.tsv", sep="\t", index=False)
T = pd.read_csv(OUT + "tables/17_protein_structure_comparison.csv")
T = T[~T.метод.astype(str).str.startswith("Boltz")]
pd.concat([T, pd.DataFrame(rows)], ignore_index=True).to_csv(OUT + "tables/17_protein_structure_comparison.csv", index=False, encoding="utf-8-sig")
