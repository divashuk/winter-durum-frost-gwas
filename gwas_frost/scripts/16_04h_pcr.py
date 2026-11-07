# -*- coding: utf-8 -*-
"""Раздел 16, шаг 4з. In silico ПЦР праймеров Pearce et al. 2013 (Table S1; tables/11_pearce2013_primers.csv) по сборкам.
Вход: blastn-short (word_size 7, evalue 1000) праймеров по синтенным окнам каждой сборки (s16/pcr_win, 16_04o_pcr_windows.sh) —
основной; по хромосомам группы 5 целиком (16_04m; Chili, Mahmoudi — вся сборка; Svevo — 5A, 5B) — дополнительный столбец.
Правила:
  * место посадки праймера: выравнивание доходит до 3′-конца праймера (qend = длина), без гэпов;
    несовпадений всего (с учётом невыровненного 5′-участка) ≤ 2; на последних 3 нуклеотидах 3′-конца несовпадений нет;
  * продукт: F и R одной пары на одной последовательности, навстречу друг другу (F на + и R на − либо наоборот),
    размер от 50 п.н. до 3000 п.н.;
  * вывод по делеции CBF-B12/B14/B15 (Fr-B2): засчитываются продукты размера ожидаемый ± 15 % в ± 2 Мб от копии CBF-B9 этой
    сборки (16_04k; так учитываются скаффолды Chili/Mahmoudi и перенос района на 7B у части мягкой пшеницы);
    «дикий тип» — есть продукты CBF-B12, CBF-B14 и CBF-B15, «делеция» — нет ни одного при наличии контроля CBF-B9;
    иначе — «частично» (перечисляются отсутствующие гены). У Chili и Mahmoudi (RagTag; скаффолды подписаны хромосомами ненадёжно)
    продукт B-специфичной пары засчитывается на любом скаффолде.
Выход: tables/16_cbf_primers_insilico.csv (пара × сборка) и s16/pcr/sites_<genome>.tsv."""
import sys, os, re
import pandas as pd, numpy as np

sys.stdout.reconfigure(encoding="utf-8")
S16 = "/work/scratch/s16/"
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qseq", "sseq"]


def sites(path):
    d = pd.read_csv(path, sep="\t", header=None, names=COLS)
    d = d[(d.qend == d.qlen) & (d.gapopen == 0)].copy()
    d["mm_total"] = d.mismatch + d.qstart - 1
    d = d[d.mm_total <= 2].copy()
    d["mm3"] = [sum(a != b for a, b in zip(q[-3:].upper(), s[-3:].upper())) for q, s in zip(d.qseq, d.sseq)]
    d = d[d.mm3 == 0].copy()
    d["pair"], d["dir"] = d.qseqid.str.split("|").str[0], d.qseqid.str.split("|").str[1]
    d["strand"] = np.where(d.send > d.sstart, "+", "-")
    d["p3"] = d.send        # координата 3′-конца праймера на мишени
    return d


def products(s, maxlen=3000):
    out = []
    for (pair, seq), g in s.groupby(["pair", "sseqid"]):
        F, R = g[g.dir == "F"], g[g.dir == "R"]
        for a, b in [(F, R), (R, F)]:      # один праймер на +, другой на −
            for _, x in a[a.strand == "+"].iterrows():
                for _, y in b[b.strand == "-"].iterrows():
                    lo, hi = x.sstart, y.sstart      # 5′-концы: x на +, y на −
                    size = hi - lo + 1
                    if 50 <= size <= maxlen:
                        out.append(dict(pair=pair, seq=seq, start=int(lo), end=int(hi), size=int(size),
                                        mm=f"{x.dir}:{x.mm_total}/{y.dir}:{y.mm_total}"))
    return pd.DataFrame(out).drop_duplicates() if out else pd.DataFrame(columns=["pair", "seq", "start", "end", "size", "mm"])


SCAF = {"Chili", "Mahmoudi"}     # RagTag: скаффолды подписаны хромосомами ненадёжно — продукт B-специфичной пары засчитывается на любом скаффолде


def seq_chr(seq):
    m = re.search(r"([1-7][ABD])", seq.replace("svevo_", ""))
    return m.group(1) if m else seq


if __name__ == "__main__":
    pr = pd.read_csv("tables/11_pearce2013_primers.csv", encoding="utf-8-sig")
    pr["pair"] = pr.gene.str.replace("-", "_").str.replace(" ", "_")
    # аргументы: <метка>=<id сборки> ; окна — s16/pcr_win/<id>.tsv, хромосомы группы 5 целиком — s16/pcr_raw/<id>.tsv (Svevo: svevo/pcr.tsv)
    files = dict(x.split("=") for x in sys.argv[1:])
    os.makedirs(S16 + "pcr", exist_ok=True)
    rows, verdict = [], []
    for g, gid in files.items():
        sw = sites(S16 + f"pcr_win/{gid}.tsv")
        Pw = products(sw)
        if len(Pw):        # окно «chr:начало-конец» -> абсолютные координаты
            m = Pw.seq.str.extract(r"^(.+):(\d+)-(\d+)$")
            off = m[1].astype(int) - 1
            Pw["chr"], Pw["start"], Pw["end"] = m[0], Pw.start + off, Pw.end + off
        else:
            Pw["chr"] = []
        fr = S16 + ("svevo/pcr.tsv" if gid == "Svevo" else f"pcr_raw/{gid}.tsv")
        Pg = products(sites(fr)) if os.path.exists(fr) and os.path.getsize(fr) else pd.DataFrame(columns=["pair", "seq", "start", "end", "size", "mm"])
        # копия CBF-B9 этой сборки (16_04k): продукты целевых пар засчитываются в ± 2 Мб от неё
        H = pd.read_csv(S16 + f"copies/{gid}_hits.tsv", sep="	") if os.path.exists(S16 + f"copies/{gid}_hits.tsv") else pd.DataFrame()
        b9 = H[(H.cand_id == "CBF_5B_523.134") & (H.ident >= 0.95)] if len(H) else H
        b9c = (b9.iloc[0].tchr, int(b9.iloc[0].tg_s)) if len(b9) else None
        for _, r in pr.iterrows():
            x = Pw[Pw.pair == r.pair].copy()
            x["size_ok"] = (x["size"] - r.expected_product_bp).abs() <= 0.15 * r.expected_product_bp
            x["near_b9"] = [b9c is not None and c == b9c[0] and abs(a - b9c[1]) <= 2e6 for c, a in zip(x.chr, x.start)]
            y = Pg[Pg.pair == r.pair]
            rows.append(dict(primer_pair=r.gene, genome=g, expected_bp=r.expected_product_bp, product="да" if len(x) else "нет",
                             n_products=len(x), sizes_bp="; ".join(map(str, x["size"])),
                             coords="; ".join(f"{c}:{a}-{b}" for c, a, b in zip(x.chr, x.start, x.end)), mismatches_F_R="; ".join(x.mm),
                             product_in_FrB2_region="да" if (x.size_ok & x.near_b9).any() else "нет",
                             products_group5_whole_chr="; ".join(f"{c}:{a}({z} п.н.)" for c, a, z in zip(y.seq, y.start, y["size"])) if len(y) else "нет"))
        R = pd.DataFrame([r for r in rows if r["genome"] == g])
        inreg = lambda p: (R[R.primer_pair.str.startswith(p)].product_in_FrB2_region == "да").any()
        tgt = {k: inreg(k) for k in ["CBF-B12", "CBF-B14", "CBF-B15"]}
        ctl = inreg("CBF-B9")
        reg = f"{b9c[0]}:{b9c[1]}" if b9c else "—"
        if all(tgt.values()) and ctl: v = "дикий тип (продукты CBF-B12, B14, B15 и контроля CBF-B9)"
        elif not any(tgt.values()) and ctl: v = "делеция (нет продуктов CBF-B12, B14, B15; контроль CBF-B9 есть)"
        elif not ctl: v = "не определено: нет продукта контроля CBF-B9"
        else: v = "частично: нет " + ", ".join(k for k, t in tgt.items() if not t)
        verdict.append(dict(genome=g, verdict=v, region_CBF_B9=reg, **{k: ("да" if t else "нет") for k, t in tgt.items()}, CBF_B9="да" if ctl else "нет"))
    out = pd.DataFrame(rows)
    V = pd.DataFrame(verdict)
    out = out.merge(V[["genome", "verdict", "region_CBF_B9"]], on="genome")
    out.to_csv("tables/16_cbf_primers_insilico.csv", index=False, encoding="utf-8-sig")
    print(V.to_string())
    print(out.pivot_table(index="primer_pair", columns="genome", values="sizes_bp", aggfunc="first").to_string())
