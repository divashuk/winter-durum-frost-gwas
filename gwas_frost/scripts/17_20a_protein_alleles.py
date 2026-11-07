# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 1а. Белковые аллели генов-кандидатов разделов 16 и 17 по 16 сборкам.
Копии — s16/copies/all_copies.tsv и s17/copies/all_copies.tsv (16_04k_integrate.py): синтенные копии своего субгенома
(идентичность ДНК ≥ 0.95, гомеолог = субгеном гена Langdon) с «ORF полная» (miniprot: покрытие белка ≥ 0.90, без сдвигов рамки и стопов).
Гены: с полной ORF хотя бы у двух сборок. Аллель — уникальная последовательность белка; A1 — аллель Langdon (если ORF Langdon полная),
далее по числу носителей. «Есть/нет»: сборки без копии в синтенном окне и сборки с копией, но неполной ORF (tables/16|17_gene_copies_by_genome.csv).
Различия пары аллелей — глобальное выравнивание белков (BLOSUM62, гэп −10/−0.5, как pair_map в 16_06_variants.py); позиция переводится
на белок Langdon; домен и консервативность позиции — из tables/16|17_variants_functional.csv (строки уровня «белок» в той же позиции p.N);
класс — правило 16_06: замена в сигнатуре CBF — сильный; в домене или в консервативной позиции при Grantham ≥ 60 — умеренный; иначе слабый;
белковый индел — умеренный (в сигнатуре — сильный).
Выход: WSL s17/prot/alleles.tsv, alleles.fa, pairs.tsv, diffs.tsv, presence.tsv, esm_jobs.tsv (позиции для ESM-1v)."""
import os, re, sys, importlib, itertools
import numpy as np, pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
v16 = lib.run16("16_06_variants.py", name="v")          # таблица Grantham и метки сборок раздела 16
GR, LABEL = v16["GRANTHAM"], v16["LABEL"]
P = lib.S17 + "prot/"; os.makedirs(P, exist_ok=True)
SVR = "участок ±20 кб вокруг гена Svevo (синтения не проверялась)"
AL = Align.PairwiseAligner(); AL.substitution_matrix = substitution_matrices.load("BLOSUM62")
AL.open_gap_score, AL.extend_gap_score, AL.mode = -10, -0.5, "global"
RANK = {"сильный": 3, "умеренный": 2, "слабый": 1}


def aligned_pairs(a, b):
    """колонки глобального выравнивания: список (i, j) с None для гэпа (1-based)"""
    aln = AL.align(a, b)[0]
    out, ia, ib = [], 0, 0
    (sa, sb) = aln[0], aln[1]
    for x, y in zip(sa, sb):
        i = (ia := ia + 1) if x != "-" else None
        j = (ib := ib + 1) if y != "-" else None
        out.append((i, j, x, y))
    return out


rows_al, rows_pair, rows_diff, rows_pres, jobs, fa = [], [], [], [], [], []
for sec, S, TC, TV, NM in [("16", lib.S16, "16_gene_copies_by_genome.csv", "16_variants_functional.csv", "candidate_names.tsv"),
                           ("17", lib.S17, "17_gene_copies_by_genome.csv", "17_variants_functional.csv", "candidate_names.tsv")]:
    if not os.path.exists(S + "copies/all_copies.tsv"): print("нет копий раздела", sec); continue
    A = pd.read_csv(S + "copies/all_copies.tsv", sep="\t")
    cand =pd.read_csv(S + "q/candidates.tsv", sep="\t").set_index("cand_id")
    names = pd.read_csv(S + "q/" + NM, sep="\t").set_index("cand_id")
    TCt = pd.read_csv(os.path.join(lib.OUT, "tables", TC)); VAR = pd.read_csv(os.path.join(lib.OUT, "tables", TV))
    for cid, x in A.groupby("cand_id", sort=False):
        if cid not in cand.index or not bool(cand.loc[cid, "is_gene"]): continue
        if sec == "17" and cid == "TRITD1Bv1G205280": continue
        sub = cand.loc[cid, "chr"][-1]
        own = x[(x.dna_ident >= 0.95) & (x.homeolog == sub) & x.synteny.isin(["в синтенном окне", SVR])].copy()
        full = own[(own.orf == "полная") & own.prot.fillna("").str.len().gt(0)].copy()
        gene = names.loc[cid, "ген"]
        lab = lambda g: LABEL.get(g, g)
        tc = TCt[TCt.cand_id == cid]
        absent = tc[tc.копий_в_синтенном_окне == 0].сборка.tolist()
        incomplete = sorted({lab(g) for g in own.genome} - {lab(g) for g in full.genome})
        rows_pres.append(dict(раздел=sec, cand_id=cid, ген=gene, сборок_с_полной_ORF=full.genome.nunique(),
                              нет_копии_в_синтенном_окне=", ".join(absent), копия_есть_ORF_неполная=", ".join(incomplete),
                              причины_неполной_ORF="; ".join(f"{lab(r.genome)}: {r.orf_why}" for _, r in own[own.orf != "полная"].iterrows())))
        if full.genome.nunique() < 2: continue
        full["p"] = full.prot.str.rstrip("*")
        full["lab"] = [lab(g) + (f"/к{k}" if k > 1 else "") for g, k in zip(full.genome, full.groupby("genome").cumcount() + 1)]
        hap = full.groupby("p").lab.apply(list)
        lang = full[full.genome == "Langdon"].p
        order = sorted(hap.index, key=lambda s: (0 if len(lang) and s == lang.iloc[0] else 1, -len(hap[s])))
        ids = {s: f"A{i + 1}" for i, s in enumerate(order)}
        Lp = x[(x.genome == "Langdon") & (x.homeolog == sub)].prot.dropna()
        Lp = Lp.iloc[0].rstrip("*") if len(Lp) else order[0]
        vp = VAR[(VAR.cand_id == cid) & (VAR.уровень == "белок")]
        posinfo = {}
        for _, r in vp.iterrows():
            m = re.match(r"p\.(\d+)", str(r.позиция_Langdon))
            if m: posinfo.setdefault(int(m.group(1)), (r.домен if isinstance(r.домен, str) else "", r.консервативность if isinstance(r.консервативность, str) else ""))
        tolang = {}
        for s in order:
            tolang[s] = {i: j for i, j, _, _ in aligned_pairs(s, Lp) if i and j}
        for s in order:
            rows_al.append(dict(раздел=sec, cand_id=cid, ген=gene, аллель=ids[s], длина=len(s), носители=", ".join(hap[s]), n_носителей=len(hap[s]),
                                как_Langdon="да" if len(lang) and s == lang.iloc[0] else "нет", белок=s))
            fa.append(f">{cid}|{ids[s]}\n{s}\n")
        for a, b in itertools.combinations(order, 2):
            cols = aligned_pairs(a, b); diffs = []; k = 0
            while k < len(cols):
                i, j, xa, yb = cols[k]
                if i and j and xa != yb:
                    L = tolang[a].get(i); dom, cons = posinfo.get(L, ("", "нет данных")) if L else ("", "нет данных")
                    g = GR.get((xa, yb))
                    cls = "сильный" if "сигнатура" in dom else ("умеренный" if dom or (cons == "консервативна" and (g or 0) >= 60) else "слабый")
                    diffs.append(dict(тип="замена", pos_a=i, aa_a=xa, pos_b=j, aa_b=yb, позиция_Langdon=L, домен=dom, консервативность=cons, Grantham=g, класс=cls))
                    k += 1
                elif i is None or j is None:
                    k0 = k; t = "вставка в " + ids[b] if i is None else "делеция в " + ids[b]
                    while k < len(cols) and ((cols[k][0] is None) == (i is None)) and ((cols[k][1] is None) == (j is None)): k += 1
                    seg = [c for c in cols[k0:k]]
                    pa = next((c[0] for c in cols[k0:] if c[0]), None)
                    L = tolang[a].get(pa) if pa else None
                    doms = {posinfo.get(tolang[a].get(c[0]), ("", ""))[0] for c in seg if c[0]} | ({posinfo.get(L, ("", ""))[0]} if L else set())
                    dom = "; ".join(sorted(d for d in doms if d))
                    diffs.append(dict(тип=f"индел ({t})", pos_a=pa, aa_a="".join(c[2] for c in seg if c[0]), pos_b=next((c[1] for c in cols[k0:] if c[1]), None),
                                      aa_b="".join(c[3] for c in seg if c[1]), позиция_Langdon=L, домен=dom, консервативность="", Grantham=None,
                                      длина=len(seg), класс="сильный" if "сигнатура" in dom else "умеренный"))
                else:
                    k += 1
            sub_ = [d for d in diffs if d["тип"] == "замена"]
            ind = [d for d in diffs if d["тип"] != "замена"]
            mx = max([RANK[d["класс"]] for d in diffs], default=0)
            rows_pair.append(dict(раздел=sec, cand_id=cid, ген=gene, a=ids[a], b=ids[b], n_замен=len(sub_), n_инделов=len(ind),
                                  макс_индел_аа=max([d["длина"] for d in ind], default=0), макс_класс={3: "сильный", 2: "умеренный", 1: "слабый", 0: "—"}[mx],
                                  макс_класс_ранг=mx))
            for d in diffs:
                rows_diff.append(dict(раздел=sec, cand_id=cid, ген=gene, a=ids[a], b=ids[b], **d))
            for d in sub_:
                jobs.append(dict(cand_id=cid, a=ids[a], b=ids[b], dir="a>b", seq_id=ids[a], pos=d["pos_a"], wt=d["aa_a"], mt=d["aa_b"]))
                jobs.append(dict(cand_id=cid, a=ids[a], b=ids[b], dir="b>a", seq_id=ids[b], pos=d["pos_b"], wt=d["aa_b"], mt=d["aa_a"]))
pd.DataFrame(rows_al).to_csv(P + "alleles.tsv", sep="\t", index=False, lineterminator="\n")
pd.DataFrame(rows_pair).to_csv(P + "pairs.tsv", sep="\t", index=False, lineterminator="\n")
pd.DataFrame(rows_diff).to_csv(P + "diffs.tsv", sep="\t", index=False, lineterminator="\n")
pd.DataFrame(rows_pres).to_csv(P + "presence.tsv", sep="\t", index=False, lineterminator="\n")
pd.DataFrame(jobs).drop_duplicates(["cand_id", "seq_id", "pos", "mt"]).to_csv(P + "esm_jobs.tsv", sep="\t", index=False, lineterminator="\n")
open(P + "alleles.fa", "w", newline="\n").write("".join(fa))
Al = pd.DataFrame(rows_al)
print("генов с ≥ 2 сборками с полной ORF:", Al.cand_id.nunique(), "аллелей:", len(Al), "пар:", len(rows_pair), "заданий ESM:", len(jobs))
print(Al.groupby("cand_id").size().describe())
