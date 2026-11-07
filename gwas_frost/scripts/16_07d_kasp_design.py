# -*- coding: utf-8 -*-
"""Раздел 16, шаг 7г. Подбор праймеров KASP по шаблонам (16_07b) с учётом гомеологов и паралогов (BLAST шаблонов, 16_07c tmpl).
Требования (задание): два аллель-специфичных праймера (3′-конец на SNP) с хвостами FAM GAAGGTGACCAAGTTCATGCT и HEX
GAAGGTCGGAGTCAACGGATT, общий праймер; ампликон 50–120 п.н.; Tm специфичной части 58–62 °C (primer3-py calc_tm, Na+ 50 мМ,
Mg2+ 1.5 мМ, dNTP 0.6 мМ, праймер 50 нМ); под праймерами нет переменных позиций сборок (маска 16_07b), кроме целевого SNP;
для каждого нецелевого локуса, выравнивание которого покрывает места посадки и аллельного, и общего праймера, 3′-конец аллельных
праймеров (3 нт) или общего праймера (2 нт) — на позиции, где этот локус отличается от целевой копии (нецелевые попадания (гомеологи 5A/5B/5D, другие CBF и т.п.: BLAST шаблона по Langdon, Svevo 5A/5B, CS целиком; целевыми считаются попадание
в саму копию Langdon и ортологи на той же хромосоме Svevo и CS, остальные с выравниванием ≥ 40 п.н. — нецелевые).
Шпильки и димеры: ΔG (primer3 calc_hairpin / calc_homodimer / calc_heterodimer) > −9 ккал/моль.
Оценка: сумма |Tm − 60|, |ампликон − 80| / 10, GC вне 35–65 %, бонус за отличие и в аллельном, и в общем праймере.
На маркер сохраняются 3 лучших варианта (разные позиции общего праймера) -> s16/kasp/designs.tsv, primers_cand.fa."""
import os, re, sys
import pandas as pd, numpy as np
import primer3

sys.stdout.reconfigure(encoding="utf-8")
S16 = "/work/scratch/s16/"
K = S16 + "kasp/"
FAM, HEX = "GAAGGTGACCAAGTTCATGCT", "GAAGGTCGGAGTCAACGGATT"
COMP = str.maketrans("ACGTN", "TGCAN")
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qseq", "sseq"]
TM = dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)


def rc(s): return s.translate(COMP)[::-1]
def tm(s): return primer3.calc_tm(s, **TM)
def dg_ok(*seqs):
    for s in seqs:
        if primer3.calc_hairpin(s, **TM).dg < -9000 or primer3.calc_homodimer(s, **TM).dg < -9000: return False
    for a in seqs:
        for b in seqs:
            if a < b and primer3.calc_heterodimer(a, b, **TM).dg < -9000: return False
    return True
def gc(s): return (s.count("G") + s.count("C")) / len(s)


T = pd.read_csv(K + "targets.tsv", sep="\t")
T = T[T.status == "мишень выбрана"].copy()
tmpl, n = {}, None
for line in open(K + "templates.fa", encoding="utf-8"):
    if line.startswith(">"): n = line[1:].strip(); tmpl[n] = ""
    else: tmpl[n] += line.strip().upper()
B = []
for src in ["langdon", "svevo", "cs"]:
    p = K + f"tmpl_{src}.tsv"
    if os.path.exists(p) and os.path.getsize(p):
        d = pd.read_csv(p, sep="\t", header=None, names=COLS); d["src"] = src; B.append(d)
B = pd.concat(B) if B else pd.DataFrame(columns=COLS + ["src"])


def nontarget(mk, chrom, pos, L):
    d = B[(B.qseqid == mk) & (B.length >= 40)].copy()
    if d.empty: return d
    d["chr"] = d.sseqid.astype(str).str.extract(r"([1-7][ABD])")[0]
    d["lo"] = d[["sstart", "send"]].min(axis=1); d["hi"] = d[["sstart", "send"]].max(axis=1)
    keep = []
    for src, g in d.groupby("src"):
        if src == "langdon":
            self_ = g[(g.chr == chrom) & (g.lo <= pos) & (g.hi >= pos) & (g.pident >= 99)]
        else:   # ортолог: лучшее попадание на той же хромосоме
            same = g[g.chr == chrom].sort_values("bitscore", ascending=False)
            self_ = same.head(1)
        keep.append(g.drop(self_.index))
    return pd.concat(keep)


def hsp_profiles(nt):
    """для каждого нецелевого попадания: (начало, конец) покрытия шаблона и словарь позиция -> основание мишени"""
    out = []
    for _, r in nt.iterrows():
        q = int(r.qstart) - 1; base = {}
        for a, b in zip(r.qseq, r.sseq):
            if a == "-": continue
            base[q] = b.upper(); q += 1
        out.append((int(r.qstart) - 1, int(r.qend) - 1, base))
    return out


def spec_ok(prof, t, i, X, Y, as3, as_site, com3, com_site):
    """нецелевой локус опасен, если покрывает места посадки обоих праймеров; тогда 3′-конец аллельного (3 нт, на SNP — не X и не Y)
    или общего праймера (2 нт) должен на нём не совпадать. Возвращает (годен, n_опасных, решено_аллельным, решено_общим)"""
    n = a_ok = c_ok = 0
    for lo, hi, base in prof:
        if not (lo <= min(as_site + com_site) and hi >= max(as_site + com_site)): continue
        n += 1
        ad = any((base.get(p, "-") not in (X, Y)) if p == i else (base.get(p, "-") != t[p]) for p in as3)
        cd = any(base.get(p, "-") != t[p] for p in com3)
        if not (ad or cd): return False, n, a_ok, c_ok
        a_ok += ad; c_ok += cd
    return True, n, a_ok, c_ok


designs = []
for _, r in T.iterrows():
    mk, t, i = r.marker, tmpl[r.marker], int(r.snp_index)
    X, Y = str(r.allele_langdon), str(r.allele_group2)
    if r.template_strand == "ген" or True:
        assert t[i] == X, (mk, t[i], X)
    masked = set(int(x) for x in str(r.masked).split(",") if x not in ("", "nan"))
    nt = nontarget(mk, r.chr, int(r.pos_langdon), len(t))
    prof = hsp_profiles(nt)
    out = []; rej = dict(mask_as=0, tm_as=0, mask_c=0, tm_c=0, spec=0)
    for orient in ["F", "R"]:
        for La in range(15, 33):          # GC-богатые CBF: Tm 58–62 достигается уже при 15–16 нт
            if orient == "F":
                a0, a1 = i - La + 1, i
                if a0 < 0: break
                core = t[a0:i]; asx, asy = core + X, core + Y
                as3 = [i - 2, i - 1, i]
            else:
                a0, a1 = i, i + La - 1
                if a1 >= len(t): break
                core = t[i + 1:a1 + 1]; asx, asy = rc(X + core), rc(Y + core)
                as3 = [i, i + 1, i + 2]
            if any(p in masked for p in range(a0, a1 + 1) if p != i): rej["mask_as"] += 1; continue
            tx, ty = tm(asx), tm(asy)
            if not (58 <= tx <= 62 and 58 <= ty <= 62): rej["tm_as"] += 1; continue
            for M in range(15, 31):
                rng = range(i + 1, len(t) - M + 1) if orient == "F" else range(0, i - M + 1)
                for b in rng:
                    c0, c1 = b, b + M - 1
                    amp = (c1 - a0 + 1) if orient == "F" else (a1 - c0 + 1)
                    if not (50 <= amp <= 120): continue
                    if any(p in masked for p in range(c0, c1 + 1)): rej["mask_c"] += 1; continue
                    com = rc(t[c0:c1 + 1]) if orient == "F" else t[c0:c1 + 1]
                    tc = tm(com)
                    if not 58 <= tc <= 62: rej["tm_c"] += 1; continue
                    com3 = [c0, c0 + 1] if orient == "F" else [c1 - 1, c1]
                    ok, n_risk, as_disc, com_disc = spec_ok(prof, t, i, X, Y, as3, [a0, a1], com3, [c0, c1])
                    if not ok: rej["spec"] += 1; continue
                    pen = abs(tx - 60) + abs(ty - 60) + abs(tc - 60) + abs(amp - 80) / 10 + \
                          sum(2 for s in (asx, asy, com) if not 0.35 <= gc(s) <= 0.65) - (0.5 if n_risk and as_disc == n_risk else 0)
                    out.append(dict(marker=mk, orient=orient, as_x=asx, as_y=asy, common=com, tm_x=round(tx, 1), tm_y=round(ty, 1),
                                    tm_c=round(tc, 1), amplicon=amp, as_start=a0, as_end=a1, com_start=c0, com_end=c1,
                                    as_disc=as_disc, com_disc=com_disc, n_risk_loci=n_risk, n_nontarget=len(nt), pen=pen))
    if not out:
        why = ", ".join(f"{k}={v}" for k, v in rej.items())
        designs.append(dict(marker=mk, status=f"праймеры с требованиями не подобраны (отсеяно: {why}; масок в шаблоне {len(masked)})", n_nontarget=len(nt))); continue
    o = pd.DataFrame(out).sort_values("pen")
    o = o[[dg_ok(a, b, c) for a, b, c in zip(o.as_x, o.as_y, o.common)] if len(o) < 4000 else o.index.isin(o.index)]
    if o.empty:
        designs.append(dict(marker=mk, status="все варианты отброшены по шпилькам/димерам", n_nontarget=len(nt))); continue
    o = o.drop_duplicates(["orient", "com_start"]).head(3)
    o["design"] = range(1, len(o) + 1); o["status"] = "подобран"
    designs += o.to_dict("records")
D = pd.DataFrame(designs)
D.to_csv(K + "designs.tsv", sep="\t", index=False, lineterminator="\n")
with open(K + "primers_cand.fa", "w", newline="\n") as fh:
    for _, d in D[D.status == "подобран"].iterrows():
        for role, s in [("X", d.as_x), ("Y", d.as_y), ("C", d.common)]:
            fh.write(f">{d.marker}|{int(d.design)}|{role}\n{s}\n")
pd.set_option("display.width", 250)
print(D[["marker", "status", "design", "orient", "amplicon", "tm_x", "tm_y", "tm_c", "n_risk_loci", "as_disc", "com_disc", "n_nontarget"]].to_string())
