# -*- coding: utf-8 -*-
"""Раздел 16, шаг 4е. Ортологи и копии генов-кандидатов по сборкам (minimap2 -x asm20, PAF шага 4в/4д).
Правила (идентичность = 1 − de, т.е. без учёта длинных инделов; инделы выписываются из CIGAR отдельно):
  * якорь найден, если лучшее попадание: идентичность ≥ 0.90, покрытие запроса ≥ 0.50;
  * синтенная хромосома группы — та, на которой больше всего найденных якорей (при равенстве — по сумме совпадений);
    синтенное окно — от крайнего до крайнего якоря этой хромосомы (якоря дальше 5 Мб от медианы отбрасываются) ± 1 Мб;
  * мишень гена — попадание запроса (ген + 2 кб + 1 кб) с идентичностью ≥ 0.90 и покрытием тела гена ≥ 0.30;
    попадания одного запроса ближе 1 кб на одной цепи сливаются; каждая мишень (по положению тела гена ± 300 п.н.)
    отдаётся запросу, в синтенное окно группы которого она попадает, а среди них — с наибольшим произведением
    идентичность × покрытие тела гена (округл. до 0.001; так разводятся паралоги кластера CBF и одинаковые гены из разных мест Langdon);
  * копия гена = назначенная мишень с идентичностью ≥ 0.95; «полная» — покрытие тела гена ≥ 0.95, иначе «неполная»;
    копия вне синтенного окна помечается «вне синтении» и в число копий в синтенном окне не входит;
  * число копий паралога = число таких копий в синтенном окне (для мягкой — по каждому гомеологу отдельно; гомеолог
    определяется по хромосоме синтенного окна группы на A, B или D — для D окно берётся по тем же якорям).
Выход: WSL s16/copies/<genome>_hits.tsv, <genome>_windows.tsv; сводная — при вызове из 16_04g."""
import re, sys, os
import numpy as np, pandas as pd

S16 = "/work/scratch/s16/"
OFF_RE = re.compile(r"^svevo_(\d[ABD])_(\d+)$")


def norm_target(t):
    """Имя мишени -> (хромосома/скаффолд, смещение)."""
    m = OFF_RE.match(t)
    if m: return m.group(1), int(m.group(2)) - 1
    if t.startswith("svevo_"): return t.replace("svevo_", ""), 0
    m = re.match(r"^(.+):(\d+)-(\d+)$", t)          # окно второго прохода «chr:начало-конец»
    if m: return m.group(1), int(m.group(2)) - 1
    return t, 0


def cig_indels(cg, min_len=10):
    out = []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg or ""):
        if op in "ID" and int(n) >= min_len: out.append(f"{op}{n}")
    return ",".join(out)


def read_paf(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        tags = dict((t[:2], t[5:]) for t in f[12:])
        tc, off = norm_target(f[5])
        rows.append(dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], tname=f[5], tchr=tc, ts=int(f[7]) + off,
                         te=int(f[8]) + off, nmatch=int(f[9]), blen=int(f[10]), tp=tags.get("tp", ""), de=float(tags.get("de", "nan")),
                         indels=cig_indels(tags.get("cg", ""))))
    d = pd.DataFrame(rows)
    if d.empty: return d
    d["ident"] = 1 - d.de
    p = d.q.str.split("|", expand=True)
    d["kind"], d["qid"], d["group"] = p[0], p[1], p[2]
    return d


def analyse(genome, paf, cand, homeolog_chrs=None):
    d = read_paf(paf)
    cand = cand.set_index("cand_id")
    A = d[d.kind == "A"].copy(); A["qcov"] = (A.qe - A.qs) / A.qlen
    A = A[(A.ident >= 0.90) & (A.qcov >= 0.5)]
    win = []
    for grp, a in A.sort_values("nmatch", ascending=False).groupby("group"):
        # для мягкой пшеницы окно ищем на каждом гомеологе: лучшие попадания якорей по каждой хромосоме группы
        a = a.copy(); a["grpchr"] = a.tchr.str.extract(r"(\d)[ABD]$")[0]
        for tchr, b in a.groupby("tchr"):
            b = b.drop_duplicates("qid")
            if len(b) < 2: continue
            med = b.ts.median(); b = b[(b.ts - med).abs() <= 5e6]
            if len(b) < 2: continue
            win.append(dict(genome=genome, group=grp, tchr=tchr, n_anchor=len(b), sum_match=b.nmatch.sum(),
                            wstart=max(0, int(b.ts.min() - 1e6)), wend=int(b.te.max() + 1e6),
                            anchors_order="+".join(b.sort_values("ts").qid.str[-6:])))
    W = pd.DataFrame(win)
    if not W.empty:
        W["rank"] = W.groupby("group").sum_match.rank(ascending=False, method="first")
    Q = d[d.kind == "Q"].copy()
    if Q.empty: return W, pd.DataFrame()
    Q["L"] = Q.qid.map(lambda c: int(cand.loc[c, "end"] - cand.loc[c, "start"] + 1))
    Q["g_s"], Q["g_e"] = 2000, 2000 + Q.L
    Q["gov"] = (np.minimum(Q.qe, Q.g_e) - np.maximum(Q.qs, Q.g_s)).clip(lower=0)
    Q["gcov"] = Q.gov / Q.L
    Q = Q[(Q.ident >= 0.90) & (Q.gcov >= 0.30)].copy()
    # положение тела гена на мишени (приблизительно, без учёта инделов)
    lead = np.maximum(Q.g_s - Q.qs, 0); tail = np.maximum(Q.qe - Q.g_e, 0)
    Q["tg_s"] = np.where(Q.strand == "+", Q.ts + lead, Q.ts + tail)
    Q["tg_e"] = np.where(Q.strand == "+", Q.te - tail, Q.te - lead)
    Q["gmatch"] = Q.gcov.clip(upper=1) * Q.ident      # счёт назначения: идентичность × покрытие тела гена
    # мишени: кластеры по положению тела гена (± 300 п.н.), назначение запросу с max совпадениями в теле гена
    Q = Q.sort_values(["tchr", "tg_s"]).reset_index(drop=True)
    tid, cur = [], (None, -1e18, 0)
    k = 0
    for _, r in Q.iterrows():
        if r.tchr != cur[0] or r.tg_s > cur[1] + 300:
            k += 1; cur = (r.tchr, r.tg_e, k)
        else:
            cur = (r.tchr, max(cur[1], r.tg_e), k)
        tid.append(k)
    Q["target"] = tid
    # синтения: мишень внутри окна группы запроса на той же хромосоме
    def syn(r):
        if W.empty: return "окна нет"
        w = W[(W.group == cand.loc[r.qid, "group"]) & (W.tchr == r.tchr)]
        return "в синтенном окне" if ((w.wstart <= r.tg_s) & (w.wend >= r.tg_e)).any() else "вне синтении"
    Q["synteny"] = Q.apply(syn, axis=1)
    # назначение: сначала запросы, в синтенное окно группы которых попадает мишень (одинаковые гены из разных мест Langdon,
    # например CBF-A12 в кластере 534 Мб и в повторе 513 Мб, иначе делили бы мишени случайно), затем по счёту
    Q["syn_ok"] = (Q.synteny == "в синтенном окне").astype(int)
    Q["gm_r"] = Q.gmatch.round(3)
    best = Q.sort_values(["target", "syn_ok", "gm_r", "gmatch"], ascending=[True, False, False, False]).groupby("target").head(1).copy()
    best["n_queries_hit"] = Q.groupby("target").qid.nunique().reindex(best.target).values
    best["second_q"] = [", ".join(sorted(set(Q[(Q.target == t) & (Q.qid != q)].qid))[:3]) for t, q in zip(best.target, best.qid)]
    best = best[best.ident >= 0.95].copy()
    best["completeness"] = np.where(best.gcov >= 0.95, "полная", "неполная")
    best["genome"] = genome
    cols = ["genome", "qid", "tchr", "tg_s", "tg_e", "strand", "ident", "gcov", "qs", "qe", "qlen", "ts", "te", "indels", "tp",
            "synteny", "completeness", "n_queries_hit", "second_q", "tname"]
    return W, best[cols].rename(columns={"qid": "cand_id"})


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    cand = pd.read_csv(S16 + "q/candidates.tsv", sep="\t")
    os.makedirs(S16 + "copies", exist_ok=True)
    for g, paf in [x.split("=") for x in sys.argv[1:]]:
        W, H = analyse(g, paf, cand)
        W.to_csv(S16 + f"copies/{g}_windows.tsv", sep="\t", index=False)
        H.to_csv(S16 + f"copies/{g}_hits.tsv", sep="\t", index=False)
        print(g, "окон", len(W), "копий", len(H))
        pd.set_option("display.width", 250)
        print(W.to_string()); print(H.drop(columns=["tname", "tp"]).to_string())
