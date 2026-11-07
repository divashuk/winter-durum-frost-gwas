# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 1в. Выбор пары аллелей на ген -> tables/17_protein_alleles.csv.
Правило (задание 17): из всех пар уникальных белков (17_20a) — пара с наибольшим классом различий (сильный > умеренный > слабый),
при равенстве — с наибольшей суммой |ESM-1v| по всем заменам между аллелями (ESM-1v — среднее 5 моделей, 17_20b). Знак суммы — направление
от аллеля a (A1 = аллель Langdon, если его ORF полная, иначе самый частый) к b: отрицательное значение — замены b менее вероятны на фоне a.
Сумма по направлению a > b — основная; b > a приводится отдельно. Инделы ESM-1v не оцениваются (только замены).
Таблица (длинный формат, столбец «строка»): «аллель» — все уникальные белки гена с носителями; «пара» — выбранная пара и её сводка;
«замена» / «индел» — различия выбранной пары: ESM-1v (5 моделей и среднее), Grantham, консервативность (ортологи ячменя, риса,
Brachypodium; tables/16|17_variants_functional.csv), домен, класс; «есть/нет» — сборки без гена или с неполной ORF.
Также s17/prot/selected_pairs.tsv — для 3D (17_21…)."""
import os, sys, importlib
import numpy as np, pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
P = lib.S17 + "prot/"
AL = pd.read_csv(P + "alleles.tsv", sep="\t"); PR = pd.read_csv(P + "pairs.tsv", sep="\t"); DF = pd.read_csv(P + "diffs.tsv", sep="\t")
PS = pd.read_csv(P + "presence.tsv", sep="\t")
E = pd.read_csv(P + "esm_scores.tsv", sep="\t").drop_duplicates(["cand_id", "seq_id", "pos", "mt"])
ek = E.set_index(["cand_id", "seq_id", "pos", "mt"])
sub = DF[DF.тип == "замена"].copy()
def esm(c, s, p, m):
    try: return ek.loc[(c, s, int(p), m)]
    except KeyError: return None
cols5 = [f"esm1v_{i}" for i in range(1, 6)]
for side, (s_col, p_col, m_col) in {"ab": ("a", "pos_a", "aa_b"), "ba": ("b", "pos_b", "aa_a")}.items():
    vals = [esm(r.cand_id, r[s_col], r[p_col], r[m_col]) for _, r in sub.iterrows()]
    sub[f"esm1v_mean_{side}"] = [v.esm1v_mean if v is not None else np.nan for v in vals]
    if side == "ab":
        for c in cols5: sub[c] = [v[c] if v is not None else np.nan for v in vals]
agg = sub.groupby(["cand_id", "a", "b"]).agg(esm_sum_ab=("esm1v_mean_ab", "sum"), esm_sum_ba=("esm1v_mean_ba", "sum"),
                                             esm_abs=("esm1v_mean_ab", lambda s: s.abs().sum())).reset_index()
PR = PR.merge(agg, on=["cand_id", "a", "b"], how="left").fillna({"esm_sum_ab": 0, "esm_sum_ba": 0, "esm_abs": 0})
best = PR.sort_values(["cand_id", "макс_класс_ранг", "esm_abs"], ascending=[True, False, False]).groupby("cand_id").head(1)
best.to_csv(P + "selected_pairs.tsv", sep="\t", index=False, lineterminator="\n")
rows = []
for _, b in best.iterrows():
    g = AL[AL.cand_id == b.cand_id]; ps = PS[PS.cand_id == b.cand_id].iloc[0]
    for _, a in g.iterrows():
        rows.append(dict(раздел=a.раздел, cand_id=a.cand_id, ген=a.ген, строка="аллель", аллель=a.аллель, длина_аа=a.длина, носители=a.носители,
                         n_носителей=a.n_носителей, как_Langdon=a.как_Langdon, выбран="да" if a.аллель in (b.a, b.b) else "нет", белок=a.белок))
    rows.append(dict(раздел=b.раздел, cand_id=b.cand_id, ген=b.ген, строка="пара", аллель=f"{b.a} vs {b.b}", выбран="да",
                     n_замен=b.n_замен, n_инделов=b.n_инделов, макс_индел_аа=b.макс_индел_аа, класс=b.макс_класс,
                     ESM1v_сумма_a_b=round(b.esm_sum_ab, 3), ESM1v_сумма_b_a=round(b.esm_sum_ba, 3), ESM1v_сумма_модулей=round(b.esm_abs, 3),
                     пар_всего=int((PR.cand_id == b.cand_id).sum())))
    rows.append(dict(раздел=b.раздел, cand_id=b.cand_id, ген=b.ген, строка="есть/нет", выбран="",
                     носители=f"нет копии в синтенном окне: {ps.нет_копии_в_синтенном_окне or '—'}; копия есть, ORF неполная: {ps.копия_есть_ORF_неполная or '—'}",
                     белок=ps.причины_неполной_ORF if isinstance(ps.причины_неполной_ORF, str) else ""))
    d = DF[(DF.cand_id == b.cand_id) & (DF.a == b.a) & (DF.b == b.b)]
    for _, r in d.iterrows():
        rec = dict(раздел=r.раздел, cand_id=r.cand_id, ген=r.ген, строка="замена" if r.тип == "замена" else "индел", аллель=f"{b.a} vs {b.b}", выбран="да",
                   вариант=(f"{r.aa_a}{int(r.pos_a)}{r.aa_b} ({b.a}→{b.b}; {b.b}: {r.aa_b}{int(r.pos_b)})" if r.тип == "замена"
                            else f"{r.тип}: {r.aa_a or ''}/{r.aa_b or ''} у {b.a}:{'' if pd.isna(r.pos_a) else int(r.pos_a)}"),
                   позиция_a=r.pos_a, позиция_b=r.pos_b, позиция_Langdon=r.позиция_Langdon, домен=r.домен if isinstance(r.домен, str) else "",
                   консервативность=r.консервативность if isinstance(r.консервативность, str) else "", Grantham=r.Grantham, класс=r.класс,
                   длина_индела=r.get("длина"))
        if r.тип == "замена":
            s = sub[(sub.cand_id == r.cand_id) & (sub.a == r.a) & (sub.b == r.b) & (sub.pos_a == r.pos_a)].iloc[0]
            rec.update({f"ESM1v_м{i}": round(s[f"esm1v_{i}"], 3) for i in range(1, 6)})
            rec.update(ESM1v_среднее_a_b=round(s.esm1v_mean_ab, 3), ESM1v_среднее_b_a=round(s.esm1v_mean_ba, 3))
        rows.append(rec)
# гены с одним аллелем (≥ 2 сборок с полной ORF, белки одинаковы)
for cid, g in AL.groupby("cand_id"):
    if cid in set(best.cand_id): continue
    ps = PS[PS.cand_id == cid].iloc[0]
    for _, a in g.iterrows():
        rows.append(dict(раздел=a.раздел, cand_id=cid, ген=a.ген, строка="аллель", аллель=a.аллель, длина_аа=a.длина, носители=a.носители,
                         n_носителей=a.n_носителей, как_Langdon=a.как_Langdon, выбран="единственный аллель — пары нет", белок=a.белок))
    rows.append(dict(раздел=ps.раздел, cand_id=cid, ген=ps.ген, строка="есть/нет",
                     носители=f"нет копии в синтенном окне: {ps.нет_копии_в_синтенном_окне or '—'}; копия есть, ORF неполная: {ps.копия_есть_ORF_неполная or '—'}"))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(lib.OUT, "tables/17_protein_alleles.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250)
print(best[["cand_id", "ген", "a", "b", "n_замен", "n_инделов", "макс_индел_аа", "макс_класс", "esm_sum_ab", "esm_abs"]].to_string())
