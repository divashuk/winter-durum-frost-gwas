# -*- coding: utf-8 -*-
"""Раздел 16, шаг 4. Локусы CBF (DREB1) по выравниваниям miniprot (GFF) белков CBF на окна сборки.
Правила:
  * выравнивание учитывается, если идентичность ≥ 0.50, выровнено ≥ 50 а.о. и длина на геноме ≤ 5 кб (гены CBF без интронов;
    miniprot запускается с -G 3k);
  * перекрывающиеся выравнивания (одна цепь) сливаются в локус;
  * лучшее попадание локуса — максимальный счёт miniprot; «паралог» — лучший ген CS той же хромосомы (5A → TraesCS5A…,
    5B → TraesCS5B…, 5D → TraesCS5D…; если такого попадания нет — лучший ген CS любой хромосомы группы 5);
    имя по T. monococcum — лучший белок Tm (Miller et al. 2006; AAY32550–AAY32564, AAX28965);
  * копия «полная», если у лучшего попадания покрытие белка ≥ 0.90, нет сдвигов рамки (Frameshift) и преждевременных стопов (StopCodon);
    иначе «неполная» (в столбце причина).
Функция call_loci используется для всех сборок; при запуске как скрипт разбирает Langdon (WSL s16/langdon_cbf/cbf.gff)."""
import re, sys
import pandas as pd

WSL = "/work/scratch/s16/"


def read_fa_len(path):
    L, n = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            n = line[1:].split()[0]; L[n] = 0
        else:
            L[n] += len(line.strip())
    return L


def parse_gff(path, plen, chr_map=None):
    rows = []
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        tname, ts, te = a["Target"].split()
        m = re.match(r"^(.+):(\d+)-(\d+)$", f[0])
        seq, off = (m.group(1), int(m.group(2)) - 1) if m else (f[0], 0)
        chrom = chr_map(seq) if chr_map else seq
        rows.append(dict(seq=seq, chr=chrom, start=int(f[3]) + off, end=int(f[4]) + off, strand=f[6], score=float(f[5]),
                         ident=float(a["Identity"]), target=tname, tstart=int(ts), tend=int(te), tlen=plen.get(tname, None),
                         fs=int(a.get("Frameshift", 0)), stop=int(a.get("StopCodon", 0))))
    d = pd.DataFrame(rows)
    d["aa"] = d.tend - d.tstart + 1
    d["cov"] = d.aa / d.tlen
    d["span"] = d.end - d.start + 1
    return d[(d.ident >= 0.5) & (d.aa >= 50) & (d.span <= 5000)].copy()


def call_loci(d, genome):
    d = d.sort_values(["chr", "strand", "start"]).reset_index(drop=True)
    lid, cur, out = 0, None, []
    for i, r in d.iterrows():
        if cur is None or r.chr != cur[0] or r.strand != cur[1] or r.start > cur[3]:
            lid += 1; cur = [r.chr, r.strand, r.start, r.end]
        else:
            cur[3] = max(cur[3], r.end)
        out.append(lid)
    d["locus"] = out
    res = []
    for k, g in d.groupby("locus"):
        b = g.loc[g.score.idxmax()]
        chrom = b.chr
        cs = g[g.target.str.startswith("CS|")]
        same = cs[cs.target.str.contains(f"TraesCS{chrom[-2:]}")] if re.search(r"[1-7][ABD]$", str(chrom)) else cs.iloc[0:0]
        csb = (same if len(same) else cs)
        csb = csb.loc[csb.score.idxmax()] if len(csb) else None
        tm = g[g.target.str.startswith("Tm|")]
        tmb = tm.loc[tm.score.idxmax()] if len(tm) else None
        sv = g[g.target.str.startswith("SV|")]
        svb = sv.loc[sv.score.idxmax()] if len(sv) else None
        why = []
        if b["cov"] < 0.9: why.append(f"покрытие {b['cov']:.2f}")
        if b.fs: why.append(f"сдвиг рамки ×{b.fs}")
        if b.stop: why.append(f"стоп ×{b.stop}")
        res.append(dict(genome=genome, chr=chrom, start=int(g.start.min()), end=int(g.end.max()), strand=b.strand,
                        best_target=b.target, best_ident=b.ident, best_cov=round(b["cov"], 3), best_fs=b.fs, best_stop=b.stop,
                        complete="полная" if not why else "неполная", why="; ".join(why),
                        cs_paralog=csb.target.split("|")[1] if csb is not None else "", cs_ident=csb.ident if csb is not None else None,
                        cs_cov=round(csb["cov"], 3) if csb is not None else None,
                        tm_name=tmb.target.split("|")[1] if tmb is not None else "", tm_ident=tmb.ident if tmb is not None else None,
                        svevo_best=svb.target.split("|")[1] if svb is not None else "", svevo_ident=svb.ident if svb is not None else None,
                        n_aln=len(g)))
    return pd.DataFrame(res)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    plen = read_fa_len(WSL + "langdon_cbf/prot.fa")
    d = parse_gff(WSL + "langdon_cbf/cbf.gff", plen, chr_map=lambda s: s.replace("langdon_", ""))
    L = call_loci(d, "Langdon")
    L.to_csv(WSL + "langdon_cbf/langdon_cbf_loci.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(L[["chr", "start", "end", "strand", "best_target", "best_ident", "best_cov", "complete", "why", "cs_paralog", "cs_ident",
             "tm_name", "tm_ident", "svevo_best", "svevo_ident"]].to_string())
