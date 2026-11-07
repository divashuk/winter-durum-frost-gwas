# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 3б (WSL, окружение s17p: TM-align, mkdssp, Biopython). Сравнение 3D-моделей двух аллелей каждой пары
(ESMFold — s17/prot/esmfold; ColabFold — s17/prot/cf, модель ранга 1) -> tables/17_protein_structure_comparison.csv.
Для пары и метода:
  * pLDDT (B-фактор CA) — среднее по белку и по каждому домену (s17/prot/domains.tsv, 17_23a);
  * TM-align (Zhang & Skolnick 2005) a против b: TM-score (нормировка по a и по b), RMSD и длина выравнивания — по всему белку и по каждому
    домену (участки доменов вырезаются из моделей; домен — по координатам на каждом аллеле);
  * вторичная структура — mkdssp (8 состояний -> H: H/G/I, E: E/B, C: прочее); доля совпадения по позициям, сопоставленным
    выравниванием белков (как в 17_20a), и состояние в каждой позиции замены у a и b;
  * для каждой замены — расстояние (CA–CA и минимальное между тяжёлыми атомами, Å) от остатка замены до ближайшего функционального
    остатка (s17/prot/functional.tsv) в модели аллеля a и b; pLDDT позиции замены.
Строки: «пара» (весь белок), «домен», «замена». Метод: ESMFold / ColabFold."""
import os, re, subprocess, sys, tempfile
import numpy as np, pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

W = "/work/scratch/s17/prot/"
OUT = "/work/project/gwas_frost/"
E = "/opt/micromamba/envs/s17p/bin/"
SP = pd.read_csv(W + "selected_pairs.tsv", sep="\t"); DM = pd.read_csv(W + "domains.tsv", sep="\t"); FN = pd.read_csv(W + "functional.tsv", sep="\t")
DF = pd.read_csv(W + "diffs.tsv", sep="\t")
seq, n = {}, None
for line in open(W + "alleles.fa"):
    if line.startswith(">"): n = line[1:].strip(); seq[n] = ""
    else: seq[n] += line.strip()
AL = Align.PairwiseAligner(); AL.substitution_matrix = substitution_matrices.load("BLOSUM62"); AL.open_gap_score, AL.extend_gap_score, AL.mode = -10, -0.5, "global"
SS3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}


def read_pdb(path):
    res = {}
    for l in open(path):
        if l.startswith("ATOM"):
            k = int(l[22:26]); el = l[76:78].strip() or l[12:16].strip()[0]
            r = res.setdefault(k, dict(atoms=[], ca=None, b=None))
            xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            if el != "H": r["atoms"].append(xyz)
            if l[12:16].strip() == "CA": r["ca"], r["b"] = xyz, float(l[60:66])
    return res


def sub_pdb(path, a, b, out):
    with open(out, "w") as f:
        for l in open(path):
            if l.startswith("ATOM") and a <= int(l[22:26]) <= b: f.write(l)
        f.write("END\n")


def tmalign(p1, p2):
    o = subprocess.run([E + "TMalign", p1, p2], capture_output=True, text=True).stdout
    g = lambda pat: (re.search(pat, o) or [None, None])[1]
    tm = re.findall(r"TM-score= ([\d.]+)", o)
    return dict(TM_a=float(tm[0]) if tm else np.nan, TM_b=float(tm[1]) if len(tm) > 1 else np.nan,
                RMSD=float(g(r"RMSD=\s*([\d.]+)") or np.nan), aligned=int(g(r"Aligned length=\s*(\d+)") or 0))


def dssp(path):
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "m.pdb")
        open(p, "w").write("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n" + open(path).read())
        o = os.path.join(d, "m.dssp")
        r = subprocess.run([E + "mkdssp", "--output-format", "dssp", p, o], capture_output=True, text=True)
        if not os.path.exists(o):
            r = subprocess.run([E + "mkdssp", p, o], capture_output=True, text=True)
        if not os.path.exists(o): return {}
        ss, on = {}, False
        for l in open(o):
            if l.startswith("  #  RESIDUE"): on = True; continue
            if on and len(l) > 16 and l[13] != "!": ss[int(l[5:10])] = SS3.get(l[16], "C")
        return ss


def cols(a, b):
    aln = AL.align(a, b)[0]; out, ia, ib = [], 0, 0
    for x, y in zip(aln[0], aln[1]):
        i = (ia := ia + 1) if x != "-" else None; j = (ib := ib + 1) if y != "-" else None
        out.append((i, j))
    return out


def mind(r1, r2):
    return float(np.min(np.linalg.norm(np.array(r1["atoms"])[:, None] - np.array(r2["atoms"])[None], axis=2)))


def model_path(method, cid, al):
    if method == "ESMFold": p = W + f"esmfold/{cid}_{al}.pdb"
    else:
        hits = sorted(x for x in os.listdir(W + "cf/out") if x.startswith(f"{cid}_{al}_") and "rank_001" in x and x.endswith(".pdb")) if os.path.isdir(W + "cf/out") else []
        p = W + "cf/out/" + hits[0] if hits else ""
    return p if p and os.path.exists(p) else None


rows = []
for _, p in SP.iterrows():
    cid, a, b = p.cand_id, p.a, p.b; sa, sb = seq[f"{cid}|{a}"], seq[f"{cid}|{b}"]
    C = cols(sa, sb); dsel = DF[(DF.cand_id == cid) & (DF.a == a) & (DF.b == b)]
    for method in ["ESMFold", "ColabFold"]:
        pa, pb = model_path(method, cid, a), model_path(method, cid, b)
        if not (pa and pb): continue
        A, B = read_pdb(pa), read_pdb(pb); ssa, ssb = dssp(pa), dssp(pb)
        pl = lambda R, lo=1, hi=10 ** 6: float(np.mean([R[k]["b"] for k in R if lo <= k <= hi])) if any(lo <= k <= hi for k in R) else np.nan
        tm = tmalign(pa, pb)
        pairs_ss = [(ssa.get(i, "C"), ssb.get(j, "C")) for i, j in C if i and j]
        rows.append(dict(ген=p.ген, cand_id=cid, аллели=f"{a} vs {b}", метод=method, строка="пара", участок="весь белок", длина_a=len(sa), длина_b=len(sb),
                         pLDDT_a=round(pl(A), 1), pLDDT_b=round(pl(B), 1), TM_score_по_a=tm["TM_a"], TM_score_по_b=tm["TM_b"], RMSD_A=tm["RMSD"],
                         выровнено_аа=tm["aligned"], совпадение_вторичной_структуры=round(np.mean([x == y for x, y in pairs_ss]), 3) if pairs_ss else np.nan,
                         вторичная_структура_a=f"H {sum(v == 'H' for v in ssa.values())}, E {sum(v == 'E' for v in ssa.values())}",
                         вторичная_структура_b=f"H {sum(v == 'H' for v in ssb.values())}, E {sum(v == 'E' for v in ssb.values())}"))
        da, db = DM[(DM.cand_id == cid) & (DM.allele == a)], DM[(DM.cand_id == cid) & (DM.allele == b)]
        for _, d in da.iterrows():
            e = db[db.domain == d.domain]
            if e.empty: continue
            e = e.iloc[0]
            with tempfile.TemporaryDirectory() as t:
                x1, x2 = os.path.join(t, "a.pdb"), os.path.join(t, "b.pdb")
                sub_pdb(pa, d.start, d.end, x1); sub_pdb(pb, e.start, e.end, x2)
                tm = tmalign(x1, x2) if (d.end - d.start) >= 4 and (e.end - e.start) >= 4 else dict(TM_a=np.nan, TM_b=np.nan, RMSD=np.nan, aligned=0)
            rows.append(dict(ген=p.ген, cand_id=cid, аллели=f"{a} vs {b}", метод=method, строка="домен", участок=f"{d.domain}: {a} {d.start}–{d.end}, {b} {e.start}–{e.end}",
                             pLDDT_a=round(pl(A, d.start, d.end), 1), pLDDT_b=round(pl(B, e.start, e.end), 1), TM_score_по_a=tm["TM_a"], TM_score_по_b=tm["TM_b"],
                             RMSD_A=tm["RMSD"], выровнено_аа=tm["aligned"]))
        fa = FN[(FN.cand_id == cid) & (FN.allele == a)]; fb = FN[(FN.cand_id == cid) & (FN.allele == b)]
        for _, s in dsel.iterrows():
            if s.тип != "замена":
                rows.append(dict(ген=p.ген, cand_id=cid, аллели=f"{a} vs {b}", метод=method, строка="индел", участок=f"{s.тип}: {a} поз. {s.pos_a} «{s.aa_a}» / {b} «{s.aa_b}»",
                                 домен=s.домен if isinstance(s.домен, str) else ""))
                continue
            i, j = int(s.pos_a), int(s.pos_b)
            def near(R, k, F):
                if k not in R or F.empty: return np.nan, np.nan, ""
                best = (1e9, 1e9, "")
                for _, f in F.iterrows():
                    if f.pos not in R or f.pos == k: continue
                    dca = float(np.linalg.norm(R[k]["ca"] - R[f.pos]["ca"])); dm = mind(R[k], R[f.pos])
                    if dm < best[1]: best = (dca, dm, f"{f.aa}{f.pos} ({f.kind})")
                if f_self := F[F.pos == k].kind.tolist(): best = (0.0, 0.0, f"сам остаток: {f_self[0]}")
                return (round(best[0], 1), round(best[1], 1), best[2]) if best[0] < 1e9 else (np.nan, np.nan, "")
            ca_a, mn_a, fr_a = near(A, i, fa); ca_b, mn_b, fr_b = near(B, j, fb)
            rows.append(dict(ген=p.ген, cand_id=cid, аллели=f"{a} vs {b}", метод=method, строка="замена", участок=f"{s.aa_a}{i}{s.aa_b} ({b}: {j})",
                             домен=s.домен if isinstance(s.домен, str) else "", класс=s.класс, pLDDT_a=A.get(i, {}).get("b"), pLDDT_b=B.get(j, {}).get("b"),
                             вторичная_структура_a=ssa.get(i, ""), вторичная_структура_b=ssb.get(j, ""),
                             до_функц_остатка_CA_a=ca_a, до_функц_остатка_мин_a=mn_a, ближайший_функц_остаток_a=fr_a,
                             до_функц_остатка_CA_b=ca_b, до_функц_остатка_мин_b=mn_b, ближайший_функц_остаток_b=fr_b))
    print(cid, a, b, flush=True)
T = pd.DataFrame(rows)
T.to_csv(OUT + "tables/17_protein_structure_comparison.csv", index=False, encoding="utf-8-sig")
print(T[T.строка == "пара"][["ген", "аллели", "метод", "pLDDT_a", "pLDDT_b", "TM_score_по_a", "RMSD_A"]].to_string())
