# -*- coding: utf-8 -*-
"""Раздел 16, шаг 4ж. Списки участков второго прохода по каждой сборке (по PAF первого прохода, правила — 16_04f_copies.py):
  win  — синтенные окна групп (все хромосомы, где найдено ≥ 2 якорей; перекрывающиеся окна сливаются) + ± 50 кб вокруг копий вне окон;
  copy — каждая копия гена-кандидата: ген + 2 кб вверх + 1 кб вниз в ориентации гена.
Выход: s16/pass2/pass2_<genome>.tsv (kind, name, chr, start, end, strand), s16/copies/<genome>_{windows,hits}.tsv.
Запуск: python scripts/16_04g_pass2_regions.py <genome>=<paf> ..."""
import sys, os, importlib.util
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("cp", os.path.join(os.path.dirname(__file__), "16_04f_copies.py"))
cp = importlib.util.module_from_spec(spec); spec.loader.exec_module(cp)
S16 = cp.S16
cand = pd.read_csv(S16 + "q/candidates.tsv", sep="\t")
os.makedirs(S16 + "pass2", exist_ok=True); os.makedirs(S16 + "copies", exist_ok=True)


def merge(iv):
    out = []
    for c, s, e in sorted(iv):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


for arg in sys.argv[1:]:
    g, paf = arg.split("=")
    W, H = cp.analyse(g, paf, cand)
    W.to_csv(S16 + f"copies/{g}_windows.tsv", sep="\t", index=False)
    H.to_csv(S16 + f"copies/{g}_hits.tsv", sep="\t", index=False)
    iv = [(r.tchr, int(r.wstart) + 1, int(r.wend)) for _, r in W.iterrows()] if len(W) else []
    for _, r in H.iterrows():
        if r.synteny != "в синтенном окне":
            iv.append((r.tchr, max(1, int(r.tg_s) - 50000), int(r.tg_e) + 50000))
    rows = [dict(kind="win", name=f"{c}:{s}-{e}", chr=c, start=s, end=e, strand="+") for c, s, e in merge(iv)]
    for _, r in H.iterrows():
        s, e = (int(r.tg_s) - 2000, int(r.tg_e) + 1000) if r.strand == "+" else (int(r.tg_s) - 1000, int(r.tg_e) + 2000)
        rows.append(dict(kind="copy", name=f"{r.cand_id}|{g}|{r.tchr}:{int(r.tg_s)}|{r.strand}", chr=r.tchr, start=max(1, s), end=e,
                         strand=r.strand))
    P = pd.DataFrame(rows)
    P.to_csv(S16 + f"pass2/pass2_{g}.tsv", sep="\t", index=False, header=False, lineterminator="\n")
    print(g, "окон", (P.kind == "win").sum(), "копий", (P.kind == "copy").sum(), "сумма окон, Мб",
          round(sum(e - s for _, _, s, e in [(0, *x) for x in merge(iv)]) / 1e6, 1) if iv else 0)
