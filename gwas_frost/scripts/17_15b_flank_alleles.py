# -*- coding: utf-8 -*-
"""Раздел 17, шаг 7з. Аллели SNP типа (б) в 16 сборках по картированию шаблонов ± 150 п.н. (17_15_server_flank_scan.sh; Svevo — локально,
участки s17/svevo и хромосомы 5A/5B ref11): лучшее попадание шаблона в сборке (наибольшее число совпадений, покрытие шаблона ≥ 0.8,
идентичность ≥ 0.95) и основание мишени против позиции SNP (151-я позиция шаблона) по строке --cs=long; «копий» — число попаданий с
идентичностью ≥ 0.95 и покрытием ≥ 0.8. Выход: s17/kasp/flank_alleles.tsv."""
import os, re, sys, glob, importlib
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
LABEL = {"ASM5281831v1": "PI 192051", "Zavitan_WEW2.1": "Zavitan", "Dicoccum_lpTriDico3": "dicoccum lpTriDico3", "CS_v1.0": "Chinese Spring",
         "Stanley": "CDC Stanley", "Landmark": "CDC Landmark", "SY_Mattis": "SY Mattis", "Norin61": "Norin 61"}
COMP = str.maketrans("ACGT", "TGCA")


rows = []
for p in sorted(glob.glob(lib.S17 + "scan_b/*.paf")) + [lib.S17 + "svevo/flanks_b_svevo.paf"]:
    if not os.path.exists(p): continue
    g = "Svevo" if "svevo" in os.path.basename(p) else os.path.basename(p)[:-4]
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        tags = dict((t[:2], t[5:]) for t in f[12:])
        q, ql, qs, qe, st = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
        nm, bl = int(f[9]), int(f[10]); cs = tags.get("cs", "")
        ident = 1 - float(tags.get("de", "1")); cov = (qe - qs) / ql
        # проход по cs (в координатах мишени; запрос — в ориентации мишени: для «−» позиции запроса идут от qe к qs)
        target = 150 if st == "+" else ql - 1 - 150
        qpos = qs if st == "+" else ql - qe        # позиция в «ориентированном» запросе
        base = None
        for op, x in re.findall(r"([=:*+\-~])([A-Za-z0-9]+)", cs):
            if op == "=":
                L = len(x)
                if qpos <= target < qpos + L: base = x[target - qpos].upper()
                qpos += L
            elif op == "*":
                if qpos == target: base = x[0].upper()          # «*ac»: a — мишень, c — запрос
                qpos += 1
            elif op == "+":
                if qpos <= target < qpos + len(x): base = "-"
                qpos += len(x)
            elif op == "-": pass
            if base is not None: break
        if base is not None and st == "-": base = base.translate(COMP)
        rows.append(dict(genome=LABEL.get(g, g), marker=q, tchr=f[5], tstart=int(f[7]), ident=round(ident, 4), qcov=round(cov, 3), nmatch=nm, allele=base))
D = pd.DataFrame(rows)
ok = D[(D.ident >= 0.95) & (D.qcov >= 0.8)]
best = ok.sort_values("nmatch", ascending=False).groupby(["marker", "genome"]).head(1)
cnt = ok.groupby(["marker", "genome"]).size().rename("копий").reset_index()
best = best.merge(cnt, on=["marker", "genome"])
best.to_csv(lib.S17 + "kasp/flank_alleles.tsv", sep="\t", index=False, lineterminator="\n")
print(best.pivot_table(index="marker", columns="genome", values="allele", aggfunc="first").to_string())
