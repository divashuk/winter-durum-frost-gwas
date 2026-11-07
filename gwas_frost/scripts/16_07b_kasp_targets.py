# -*- coding: utf-8 -*-
"""Раздел 16, шаг 7б. Мишени KASP и шаблоны (Langdon ± 150 п.н.).
(а) ген-кандидат (приоритет 1 и 2, кроме фрагментов и повтора 5A 513 Мб): колонка геномного выравнивания (genes16/aln/<ген>_genomic.aln.fa,
    только синтенные копии), где все сборки группы 2 (16_haplotype_groups.csv) несут один нуклеотид, а все сборки группы 1 — другой,
    без гэпов в колонке; приоритет: CDS/тело гена, затем промотор 1 кб, затем фланки; при равенстве — вариант с наибольшим классом
    (16_variants_functional.csv), затем ближе к середине гена. Если группа 2 отличается только отсутствием гена — KASP не применим.
(б) SNP GWAS: 5B:535723924, 5B:535618213, 5B:522742783, 5B:532482598 (аллели — VCF коллекции).
Маскируются (под праймеры не попадают) все переменные колонки выравнивания в шаблоне, кроме целевой; для (б) — по BLAST фланков
по сборкам (16_04m_server_pcr.sh) и SNP VCF коллекции в ± 150 п.н.
Выход: s16/kasp/targets.tsv, templates.fa."""
import os, re, sys
import pandas as pd, numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
S16 = "/work/scratch/s16/"
sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location("v", os.path.join(HERE, "16_06_variants.py")); v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
COMP = str.maketrans("ACGTN", "TGCAN")
FL = 150


def rc(s): return s.translate(COMP)[::-1]


names = pd.read_csv(S16 + "q/candidate_names.tsv", sep="\t")
fb = pd.read_csv(S16 + "q/file_base.tsv", sep="\t", header=None, names=["base", "cand_id"])
H = pd.read_csv(os.path.join(OUT, "tables/16_haplotype_groups.csv"))
VAR = pd.read_csv(os.path.join(OUT, "tables/16_variants_functional.csv"))
cand = pd.read_csv(S16 + "q/candidates.tsv", sep="\t").set_index("cand_id")
rows, tmpl = [], []
for _, f in fb.iterrows():
    nm = names[names.cand_id == f.cand_id].iloc[0]
    if not str(nm.приоритет).startswith(("1", "2")) or "повтор" in str(nm.приоритет): continue
    h = H[H.cand_id == f.cand_id]
    if h.empty: continue
    h = h.iloc[0]
    g1 = set(str(h.группа_1_как_Langdon).split(", ")) if isinstance(h.группа_1_как_Langdon, str) else set()
    g2 = set(str(h.группа_2).split(", ")) if isinstance(h.группа_2, str) else set()
    gn, gs = v.read_aln(os.path.join(v.ALN, f"{f.base}_genomic.aln.fa"))
    ri = v.ref_index(gn)
    if ri is None: continue
    lab = [n.split("|")[0].replace("_", " ") for n in gn]
    in1 = [k for k, l in enumerate(lab) if l in g1]; in2 = [k for k, l in enumerate(lab) if l in g2]
    c = cand.loc[f.cand_id]; L = int(c.end - c.start + 1)
    rec = dict(marker="", type="а: ген-кандидат", gene=nm.ген, cand_id=f.cand_id, group1=", ".join(sorted(g1)), group2=", ".join(sorted(g2)),
               group2_diff=h.отличие_группы_2)
    if not in2:
        rows.append({**rec, "status": "нет сборок группы 2 с последовательностью гена (" + str(h.отличие_группы_2) + ") — KASP не применим"}); continue
    refc = v.colmap(gs[ri]); cands = []
    var_cols = [k for k in range(len(gs[ri])) if len({s[k] for s in gs if s[k] != "-"}) > 1 or any(s[k] == "-" for s in gs)]
    for k in range(len(gs[ri])):
        a1 = {gs[j][k] for j in in1 + [ri]}; a2 = {gs[j][k] for j in in2}
        if len(a1) != 1 or len(a2) != 1 or a1 == a2 or "-" in a1 | a2 or "N" in a1 | a2: continue
        near_gap = any(any(s[x] == "-" for s in gs) for x in range(max(0, k - 10), min(len(gs[ri]), k + 11)))
        p = refc[k]
        zone = 0 if 2000 < p <= 2000 + L else (1 if 1000 < p <= 2000 else 2)
        cands.append(((zone, near_gap, abs(p - (2000 + L / 2))), k, p, next(iter(a1)), next(iter(a2))))
    if not cands:
        rows.append({**rec, "status": "нет колонки, полностью разделяющей группы (SNP без гэпов)"}); continue
    picked = []
    for cnd in sorted(cands):                 # до 5 мишеней на ген, не ближе 20 п.н. друг к другу (по порядку предпочтения)
        if all(abs(cnd[2] - q[2]) >= 20 for q in picked): picked.append(cnd)
        if len(picked) == 5: break
    ref = gs[ri].replace("-", "")
    for rank, (_, k, p, x, y) in enumerate(picked, 1):
        lo, hi = max(0, p - 1 - FL), min(len(ref), p + FL)
        t = ref[lo:hi]; snp_i = p - 1 - lo
        masked = sorted({refc[q] - 1 - lo for q in var_cols if refc[q] and q != k and lo <= refc[q] - 1 < hi})
        rel = f"g.{p - 2001}" if p <= 2000 else (f"g.{p - 2000}" if p <= 2000 + L else f"g.*{p - 2000 - L}")
        # координата Langdon: копия Langdon — в ориентации гена; ген + 2 кб вверх + 1 кб вниз
        lc = int(c.qstart) + p - 1 if c.strand == "+" else int(c.qend) - p + 1
        zone = "CDS / тело гена" if 2000 < p <= 2000 + L else ("промотор 1 кб" if 1000 < p <= 2000 else "фланк")
        vv = VAR[(VAR.cand_id == f.cand_id) & (VAR.позиция_Langdon == rel)]
        alleles_by = {lab[j]: gs[j][k] for j in range(len(gs))}
        mk = f"KASP_{re.sub(r'[^A-Za-z0-9]+', '_', nm.ген)[:28]}_{c.chr}_{lc}"
        rows.append({**rec, "marker": mk, "target_rank": rank, "status": "мишень выбрана", "chr": c.chr, "pos_langdon": lc, "rel": rel, "zone": zone,
                     "allele_langdon": x, "allele_group2": y,
                     "variant_note": "; ".join(vv.вариант + " (" + vv.эффект + ", " + vv.класс + ")") if len(vv) else "",
                     "alleles_assemblies": "; ".join(f"{a}:{b}" for a, b in alleles_by.items()), "template_strand": "ген", "snp_index": snp_i,
                     "masked": ",".join(map(str, masked))})
        tmpl.append((mk, t))
# (б) SNP GWAS
pg = pd.read_csv(S16 + "kasp/parents_gt.tsv", sep="\t")
pg.columns = [re.sub(r"^#?\[\d+\]", "", c).replace(":GT", "") for c in pg.columns]
fl = {}
name = None
for line in open(S16 + "kasp/gwas_flanks.fa", encoding="utf-8"):
    if line.startswith(">"): name = line[1:].strip(); fl[name] = ""
    else: fl[name] += line.strip().upper()
for snp in ["535723924", "535618213", "522742783", "532482598"]:
    r = pg[pg.POS == int(snp)].iloc[0]
    t = fl[f"SNP_5B_{snp}"]
    assert t[FL] == r.REF, (snp, t[FL], r.REF)
    mk = f"KASP_5B_{snp}"
    rows.append(dict(marker=mk, type="б: SNP GWAS", gene="", cand_id="", target_rank=1, status="мишень выбрана", chr="5B", pos_langdon=int(snp), rel="",
                     zone="SNP GBS", allele_langdon=r.REF, allele_group2=r.ALT, variant_note="SNP GBS коллекции (FINAL-181124)",
                     template_strand="+", snp_index=FL, masked=""))
    tmpl.append((mk, t))
T = pd.DataFrame(rows)
T.to_csv(S16 + "kasp/targets.tsv", sep="\t", index=False, lineterminator="\n")
with open(S16 + "kasp/templates.fa", "w", newline="\n") as fh:
    for n, s in tmpl: fh.write(f">{n}\n{s}\n")
pd.set_option("display.width", 250)
print(T[["marker", "gene", "status", "zone", "allele_langdon", "allele_group2", "group2"]].to_string())
