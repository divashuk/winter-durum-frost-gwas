# -*- coding: utf-8 -*-
"""Этап 10.5. Пересечение локусов с опубликованными QTL/MTA морозо- и холодостойкости пшеницы.
Вход: таблица литературы (584 строки, 34 источника; собрана по статьям, DOI проверены в Crossref) —
копия сохраняется как tables/10_literature_qtl_source.csv.
Сопоставление без пересчёта координат между сборками (цепочек liftover нет), по координатной системе источника:
 * Langdon: маркеры нашей статьи и Divashuk et al. 2026 Plants — положения по последовательности праймеров/генов (этапы 5 и 9.4);
 * Svevo v1: границы локуса в Svevo v1 = положения генов Svevo, перенесённых в широкий интервал (5–95 % при >= 10 генах, иначе min–max);
 * IWGSC RefSeq v1.0 (Ensembl Plants, GCA_900519105.1): границы = положения лучших ортологов CS тех же генов на той же хромосоме;
 * RefSeq v2.1, «не указано», CerealsDB, WheatOmics, Svevo v2.0 — сравниваются с RefSeq v1.0 / Svevo v1 ПРИБЛИЖЁННО (помечено);
 * только сМ или без позиции — сопоставление по названию известного локуса (Fr-A2, Fr-B2, VRN-A1, VRN-B1, PPD-A1), иначе «не сопоставимо».
Связь: «перекрывается» (интервалы пересекаются) или «рядом» (до 5 Мб) — с указанием расстояния."""
import os, re, sys
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
W = r"/work/scratch/genes10"
T = r"/work/project/gwas_frost/tables"
SRC = r"/work/project/gwas_frost/tables/lit_frost_qtl.csv"
lit = pd.read_csv(SRC)
lit.to_csv(os.path.join(T, "10_literature_qtl_source.csv"), index=False, encoding="utf-8-sig")
L = pd.read_csv(os.path.join(T, "10_loci.csv"))
G = pd.read_csv(os.path.join(T, "10_genes_wide.csv"))


def span(gg, col):
    # только гены, коллинеарные Langdon: смещение (Langdon − другая сборка) отличается от медианного не более чем на 5 Мб
    # (паралоги/ортологи в другом месте хромосомы иначе раздувают границы)
    x = gg[["langdon_start", col]].dropna().astype(float)
    if x.empty:
        return (np.nan, np.nan)
    off = x.langdon_start - x[col]; x = np.sort(x.loc[(off - off.median()).abs() <= 5e6, col].values)
    return (np.quantile(x, .05), np.quantile(x, .95)) if len(x) >= 10 else (x[0], x[-1])


sp = []
for r in L.itertuples():
    g = G[G.locus_id == r.locus_id]
    sv = span(g[g.svevo_chr == r.chr], "svevo_start"); cs = span(g[g.cs_chr.astype(str) == r.chr], "cs_start")
    gn = g[g.in_narrow]; svn = span(gn[gn.svevo_chr == r.chr], "svevo_start"); csn = span(gn[gn.cs_chr.astype(str) == r.chr], "cs_start")
    sp.append(dict(locus_id=r.locus_id, name=r.name, chr=r.chr, langdon=(r.qtl_start, r.qtl_end), langdon_n=(r.narrow_start, r.narrow_end),
                   svevo=sv, svevo_n=svn, cs=cs, cs_n=csn))
SP = pd.DataFrame(sp)
pd.DataFrame([{"locus_id": s.locus_id, "locus_name": s.name, "chr": s.chr,
               "langdon_wide": f"{s.langdon[0] / 1e6:.2f}–{s.langdon[1] / 1e6:.2f}", "langdon_narrow": f"{s.langdon_n[0] / 1e6:.2f}–{s.langdon_n[1] / 1e6:.2f}",
               "svevo_v1_wide": f"{s.svevo[0] / 1e6:.2f}–{s.svevo[1] / 1e6:.2f}", "svevo_v1_narrow": f"{s.svevo_n[0] / 1e6:.2f}–{s.svevo_n[1] / 1e6:.2f}",
               "cs_v1_wide": f"{s.cs[0] / 1e6:.2f}–{s.cs[1] / 1e6:.2f}", "cs_v1_narrow": f"{s.cs_n[0] / 1e6:.2f}–{s.cs_n[1] / 1e6:.2f}"} for s in SP.itertuples()]
             ).to_csv(os.path.join(T, "10_loci_coordinates_svevo_cs.csv"), index=False, encoding="utf-8-sig")

LANGDON_MARKER = {"S2269949": 524288011, "S1077313": 537231599, "S1862541": 530965192, "S1298957": 513240733, "S1051014": 549361380,
                  "Mor534": 550225979, "5B_517276534": 550225979, "Mor587": 54736990, "1B_41099587": 54736990, "Mor866": 610783371, "7B_598228866": 610783371,
                  "Vrn1_new": 598268462, "VrnA1_9K0001": 598268462, "Vrn-B1": 605129460}
NAME_RULES = [(r"Fr-?B2|CBF-B\d", ["L22"]), (r"Fr-?A2|CBF-A\d|CBF14|CBF-A14", ["L02", "L01"]), (r"VRN-?A1|Vrn-?A1|Fr-A1", ["L23"]),
              (r"VRN-?B1|Vrn-?B1", ["L09"]), (r"PPD-?A1|Ppd-?A1", ["L11"])]


def rel(s, e, a, b):
    if np.isnan(a):
        return None
    if s <= b and e >= a:
        return 0.0
    return min(abs(s - b), abs(a - e)) / 1e6


rows = []
for r in lit.itertuples():
    out = dict(r._asdict()); out.pop("Index")
    s, e = r.pos_start_bp, r.pos_end_bp
    asm = str(r.assembly)
    mk = [m for m in LANGDON_MARKER if isinstance(r.marker, str) and re.search(r"(^|[^\w])" + re.escape(m) + r"($|[^\w])", r.marker)]
    if not mk and isinstance(r.locus_name, str):
        mk = [m for m in LANGDON_MARKER if m in r.locus_name and m.startswith(("1B_", "5B_", "7B_"))]
    hits, method, approx = [], "", False
    if mk:
        pos = LANGDON_MARKER[mk[0]]; s = e = pos; method = f"Langdon по последовательности маркера {mk[0]} (этапы 5/9.4)"; key = "langdon"
    elif asm.startswith("genetic cM") or (pd.isna(s)):
        key = None; method = "только сМ" if asm.startswith("genetic cM") else "позиции нет"
    elif asm.startswith("Svevo v1"):
        key = "svevo"; method = "Svevo v1 (границы локуса по перенесённым генам)"
    elif asm.startswith("Svevo v2"):
        key = "svevo"; approx = True; method = "Svevo v2.0 сравнивается с Svevo v1 приближённо"
    elif asm == "IWGSC RefSeq v1.0" or asm.startswith("IWGSC RefSeq v1.0 ("):
        key = "cs"; method = "IWGSC RefSeq v1.0 (границы локуса по ортологам CS)"
    else:
        key = "cs"; approx = True; method = f"{asm} сравнивается с RefSeq v1.0 приближённо"
    if key:
        for q in SP[SP.chr == str(r.chr)].itertuples():
            a, b = getattr(q, key); an, bn = getattr(q, key + "_n")
            d = rel(s, e, a, b)
            if d is not None and d <= 5:
                dn = rel(s, e, an, bn)
                hits.append(f"{q.locus_id} {q.name}: {'перекрывается' if d == 0 else f'рядом, {d:.1f} Мб'}"
                            + ("" if dn is None else f"; узкий: {'внутри' if dn == 0 else f'{dn:.1f} Мб'}"))
    if not hits and isinstance(r.locus_name, str):
        for pat, ids in NAME_RULES:
            if re.search(pat, r.locus_name + " " + str(r.marker)) and SP[SP.locus_id.isin(ids)].chr.eq(str(r.chr)).any():
                hits = [f"{i} {SP.loc[SP.locus_id == i, 'name'].iloc[0]}: по названию локуса" for i in ids]; method += "; сопоставлено по названию"
                break
    out.update(comparison=method, approximate=approx, matched_loci=" | ".join(hits),
               relation=("перекрывается" if any("перекрывается" in h for h in hits) else "рядом (<= 5 Мб)" if any("рядом" in h for h in hits)
                         else "по названию" if hits else "не сопоставимо" if not key else "нет пересечения"))
    rows.append(out)
O = pd.DataFrame(rows)
O.to_csv(os.path.join(T, "10_literature_overlap.csv"), index=False, encoding="utf-8-sig")
print(O.relation.value_counts().to_string())
# сводка по локусам
summ = []
for q in SP.itertuples():
    m = O[O.matched_loci.str.contains(q.locus_id + " ", regex=False)]
    summ.append(dict(locus_id=q.locus_id, locus_name=q.name, n_overlap=int(m.matched_loci.str.contains(q.locus_id + r" [^|]*перекрывается", regex=True).sum()),
                     n_near=int(m.matched_loci.str.contains(q.locus_id + r" [^|]*рядом", regex=True).sum()),
                     n_name=int(m.matched_loci.str.contains(q.locus_id + r" [^|]*по названию", regex=True).sum()),
                     sources=", ".join(sorted(m.source_id.unique()))))
S = pd.DataFrame(summ); S.to_csv(os.path.join(T, "10_literature_overlap_by_locus.csv"), index=False, encoding="utf-8-sig")
print(S.to_string())
