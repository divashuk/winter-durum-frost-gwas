# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 3а. Домены и функциональные остатки на каждом аллеле выбранных пар (s17/prot/selected_pairs.tsv).
Домены: InterPro (Ensembl Plants BioMart, координаты на белке) белка CS-паралога (CBF) или Svevo (прочие гены), перенесённые на белок аллеля
глобальным выравниванием (pair_map из 16_06_variants.py); записи длиной ≥ 80 % белка (семейство, суперсемейство) доменом не считаются;
для CBF — ещё сигнатуры PKK…ETRHP и DSAWR (find_sig из 16_06).
Функциональные остатки:
  CBF — остатки AP2, контактирующие с ДНК в комплексе AtERF1–GCC-бокс (PDB 1GCC, Allen et al. 1998; модель 1 ЯМР-ансамбля): тяжёлый атом
        белка в ≤ 4.5 Å от тяжёлого атома ДНК; перенос на аллель — локальное выравнивание белка 1GCC с аллелем (BLOSUM62);
  PRR95 — фосфоакцепторная позиция receiver-домена: остаток, выровненный с D57 CheY E. coli (UniProt P0AE67; глобальное выравнивание
        receiver-домена аллеля, InterPro IPR001789, с CheY); CCT — остатки домена IPR010402;
  ферменты: TRITD5Av1G171440 (GH1) — каталитические Glu по мотивам PROSITE PS00572 (нуклеофил, [LIVMFSTC][LIVFYS][LIV][LIVMST]ENG) и
        T[FLY]NEP (кислота/основание GH1); NB-ARC (TRITD1Bv1G205280) — P-петля (Walker A, G-x(4)-GK[ST]);
  прочие гены — функциональные остатки не определялись; расстояние даётся до ближайшего остатка домена InterPro (отмечено в таблице).
Выход: s17/prot/domains.tsv, functional.tsv; s17/prot/ref/1GCC.pdb, P0AE67.fasta."""
import os, re, sys, importlib, urllib.request
import numpy as np, pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
v = lib.run16("16_06_variants.py", name="v")
P = lib.S17 + "prot/"; R = P + "ref/"; os.makedirs(R, exist_ok=True)
for url, f in [("https://files.rcsb.org/download/1GCC.pdb", "1GCC.pdb"), ("https://rest.uniprot.org/uniprotkb/P0AE67.fasta", "P0AE67.fasta")]:
    if not os.path.exists(R + f): urllib.request.urlretrieve(url, R + f)
SP = pd.read_csv(P + "selected_pairs.tsv", sep="\t")
seq, n = {}, None
for line in open(P + "alleles.fa"):
    if line.startswith(">"): n = line[1:].strip(); seq[n] = ""
    else: seq[n] += line.strip()
three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
st = PDBParser(QUIET=True).get_structure("1gcc", R + "1GCC.pdb")[0]
prot_res = [r for ch in st for r in ch if r.get_resname() in three]
dna_atoms = np.array([a.coord for ch in st for r in ch if r.get_resname().strip() in ("DA", "DC", "DG", "DT") for a in r if a.element != "H"])
erf_seq = "".join(three[r.get_resname()] for r in prot_res)
contact = []
for k, r in enumerate(prot_res):
    xyz = np.array([a.coord for a in r if a.element != "H"])
    if len(xyz) and np.min(np.linalg.norm(xyz[:, None] - dna_atoms[None], axis=2)) <= 4.5:
        contact.append((k + 1, three[r.get_resname()] + str(r.id[1])))
print("1GCC: белок", len(erf_seq), "а.о.; контактов с ДНК", len(contact), [c[1] for c in contact])
chey = "".join(l.strip() for l in open(R + "P0AE67.fasta") if not l.startswith(">"))
LA = Align.PairwiseAligner(); LA.substitution_matrix = substitution_matrices.load("BLOSUM62"); LA.open_gap_score, LA.extend_gap_score, LA.mode = -10, -0.5, "local"


def amap(a, b, mode="global"):
    if mode == "global": return v["pair_map"](a, b)
    aln = LA.align(a, b)[0]; m = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1): m[s1 + k + 1] = s2 + k + 1
    return m


src = {}
for sec, S in [("16", lib.S16), ("17", lib.S17)]:
    nm = pd.read_csv(S + "q/candidate_names.tsv", sep="\t").set_index("cand_id")
    ipr_cs = pd.read_csv(S + "biomart/cs_interpro.tsv", sep="\t").drop_duplicates(["ensembl_gene_id", "interpro", "interpro_start"])
    ipr_sv = pd.read_csv(S + "biomart/svevo_interpro.tsv", sep="\t").drop_duplicates(["ensembl_gene_id", "interpro", "interpro_start"])
    cs_pep = {re.search(r"gene:(\S+)", h).group(1): s for h, s in v["read_fa_full"](S + "seq/cs_pep.fa").items()}
    sv_pep = {h.split()[0]: s for h, s in v["read_fa_full"](S + "seq/svevo_pep.fa").items()}
    src[sec] = (nm, ipr_cs, ipr_sv, cs_pep, sv_pep)
D, F = [], []
for _, p in SP.iterrows():
    nm, ipr_cs, ipr_sv, cs_pep, sv_pep = src[str(p.раздел)]
    r = nm.loc[p.cand_id]; is_cbf = p.cand_id.startswith("CBF_")
    sid = r.ID_CS_RefSeq_v1_1 if is_cbf else str(r.ID_Svevo).split(" ")[0]
    pep = (cs_pep if is_cbf else sv_pep).get(sid); ipr = (ipr_cs if is_cbf else ipr_sv)
    ipr = ipr[ipr.ensembl_gene_id == sid].dropna(subset=["interpro_start"])
    for al in (p.a, p.b):
        s = seq[f"{p.cand_id}|{al}"]; doms = []
        if pep:
            m = amap(pep.rstrip("*"), s)
            for _, d in ipr.iterrows():
                a_, b_ = m.get(int(d.interpro_start)), m.get(int(d.interpro_end))
                if a_ and b_ and (b_ - a_ + 1) < 0.8 * len(s): doms.append((f"{d.interpro_description} ({d.interpro})", a_, b_))
        if is_cbf:
            for k, sig in v["SIG"].items():
                x = v["find_sig"](s, sig)
                if x: doms.append((k, x[0], x[1]))
        seen = set()
        for dn, a_, b_ in doms:
            if (dn, a_, b_) in seen: continue
            seen.add((dn, a_, b_)); D.append(dict(cand_id=p.cand_id, allele=al, domain=dn, start=a_, end=b_))
        fr = []
        if is_cbf:
            m = amap(erf_seq, s, "local")
            fr = [(m[k], f"контакт с ДНК в 1GCC: AtERF1 {lab}") for k, lab in contact if k in m]
        elif p.cand_id == "TRITD5Bv1G171760":
            rec = [(a_, b_) for dn, a_, b_ in doms if "IPR001789" in dn]
            if rec:
                a_, b_ = min(x[0] for x in rec), max(x[1] for x in rec)
                m = amap(chey, s[a_ - 1:b_])
                if 57 in m: fr.append((a_ - 1 + m[57], f"позиция фосфоакцепторного D57 CheY (у аллеля {s[a_ - 2 + m[57]]})"))
            for dn, a_, b_ in doms:
                if "IPR010402" in dn: fr += [(k, "домен CCT (IPR010402)") for k in range(a_, b_ + 1)]
        elif p.cand_id == "TRITD5Av1G171440":
            for mo in re.finditer(r"(?=([LIVMFSTC][LIVFYS][LIV][LIVMST]ENG))", s): fr.append((mo.start() + 5, "каталитический нуклеофил GH1 (мотив PS00572 …ENG)"))
            for mo in re.finditer(r"(?=(T[FLY]NEP))", s): fr.append((mo.start() + 4, "кислота/основание GH1 (мотив T[FLY]NEP)"))
        elif p.cand_id == "TRITD1Bv1G205280":
            for mo in re.finditer(r"(?=(G.{4}GK[ST]))", s): fr += [(mo.start() + k + 1, "P-петля NB-ARC (Walker A)") for k in range(8)]
        if not fr:
            fr = [(k, f"остаток домена {dn} (функциональные остатки не определялись)") for dn, a_, b_ in doms if "сигнатура" not in dn for k in range(a_, b_ + 1)]
        for k, lab in fr:
            F.append(dict(cand_id=p.cand_id, allele=al, pos=int(k), aa=s[int(k) - 1], kind=lab))
pd.DataFrame(D).to_csv(P + "domains.tsv", sep="\t", index=False, lineterminator="\n")
pd.DataFrame(F).drop_duplicates(["cand_id", "allele", "pos", "kind"]).to_csv(P + "functional.tsv", sep="\t", index=False, lineterminator="\n")
print("доменов", len(D), "функциональных остатков", len(F))
