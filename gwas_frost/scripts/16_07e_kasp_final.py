# -*- coding: utf-8 -*-
"""Раздел 16, шаг 7д. Окончательный выбор KASP и таблицы tables/16_kasp_primers.csv, tables/16_kasp_specificity.csv.
Попадание праймера (blastn-short по Langdon — 14 хромосом, Svevo — 5A и 5B, CS RefSeq v1.0 — весь геном): выравнивание доходит до
3′-конца, без гэпов, несовпадений ≤ 2 (с учётом невыровненного 5′-участка), последние 3 нт совпадают.
Целевое попадание: Langdon — в ± 200 п.н. от SNP; Svevo и CS — внутри лучшего попадания шаблона на той же хромосоме (ортолог).
Посторонний продукт: аллельный и общий праймер садятся на одну последовательность навстречу друг другу не дальше 1000 п.н.
(не в целевом месте). Выбор на ген: без посторонних продуктов ни в одной из трёх сборок и с целевым продуктом в Langdon, затем
наименьшее число одиночных нецелевых посадок, затем ранг мишени (16_07b) и вариант праймеров (16_07d).
Хвосты: FAM GAAGGTGACCAAGTTCATGCT — аллель Langdon; HEX GAAGGTCGGAGTCAACGGATT — второй аллель.
Ожидаемые генотипы: Langdon и сборки — по выравниванию (16_07b) или фланку (тип б); родители F2 и SVEVO коллекции — по GBS,
только если SNP маркера — SNP VCF коллекции; иначе «нет данных GBS».
Делеция CBF-B12/B14/B15 — строка доминантной ПЦР (праймеры Pearce et al. 2013, выбор 12_pcr_primer_choice.csv) с ожидаемыми
продуктами по 16_cbf_primers_insilico.csv и GBS-гаплотипом Fr-B2 родителей (14_f2_parents_genotypes.csv, без предсказания продукта)."""
import os, re, sys
import pandas as pd, numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
S16 = "/work/scratch/s16/"; K = S16 + "kasp/"
FAM, HEX = "GAAGGTGACCAAGTTCATGCT", "GAAGGTCGGAGTCAACGGATT"
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qseq", "sseq"]
PAR = [("TP_MO_081", "Цель"), ("TP_MO_044", "Андромеда"), ("TP_MO_054", "Прибуткова"), ("TP_MO_014", "Синьора"), ("TP_MO_005", "Алена"),
       ("TP_MO_037", "Степной янтарь"), ("TP_MO_171", "ANTALIS"), ("TP_MO_190", "PAPRICA"), ("TP_MO_196", "CORTEZ"), ("TP_MO_302", "PALEOTTO"),
       ("TP_MO_304", "L2870 VILLOSA"), ("TP_MO_350", "RAMIRES")]
GEN = {"langdon": "Langdon", "svevo": "Svevo (5A, 5B)", "cs": "Chinese Spring RefSeq v1.0"}


def load(p):
    if not os.path.exists(p) or not os.path.getsize(p): return pd.DataFrame(columns=COLS)
    return pd.read_csv(p, sep="\t", header=None, names=COLS)


T = pd.read_csv(K + "targets.tsv", sep="\t")
D = pd.read_csv(K + "designs.tsv", sep="\t")
D = D[D.status == "подобран"].copy(); D["design"] = D.design.astype(int)
PH = {g: load(K + f"primer_{g}.tsv") for g in GEN}
TH = {g: load(K + f"tmpl_{g}.tsv") for g in GEN}
hits = []
for g, d in PH.items():
    if d.empty: continue
    d = d[(d.qend == d.qlen) & (d.gapopen == 0)].copy()
    d["mm"] = d.mismatch + d.qstart - 1
    d = d[d.mm <= 2]
    d = d[[q[-3:].upper() == s[-3:].upper() for q, s in zip(d.qseq, d.sseq)]]
    d["genome"] = g; d["strand"] = np.where(d.send > d.sstart, "+", "-")
    d["chr"] = d.sseqid.astype(str).str.extract(r"([1-7][ABD]|Un)")[0]
    hits.append(d)
H = pd.concat(hits) if hits else pd.DataFrame(columns=COLS + ["genome", "strand", "chr", "mm"])
H[["marker", "design", "role"]] = H.qseqid.str.split("|", expand=True) if len(H) else None


def target_span(mk, g, chrom, pos):
    if g == "langdon": return (f"langdon_{chrom}", pos - 200, pos + 200)
    t = TH[g]; t = t[t.qseqid == mk].copy()
    if t.empty: return None
    t["chr"] = t.sseqid.astype(str).str.extract(r"([1-7][ABD])")[0]
    t = t[t.chr == chrom].sort_values("bitscore", ascending=False)
    if t.empty: return None
    r = t.iloc[0]; return (r.sseqid, min(r.sstart, r.send) - 50, max(r.sstart, r.send) + 50)


rows_spec, eval_ = [], []
for _, d in D.iterrows():
    t = T[T.marker == d.marker].iloc[0]
    res = dict(marker=d.marker, design=d.design, off_products=0, off_hits=0, target_langdon=False)
    for g in GEN:
        span = target_span(d.marker, g, t.chr, int(t.pos_langdon))
        hx = H[(H.marker == d.marker) & (H.design == str(d.design)) & (H.genome == g)]
        is_t = hx.apply(lambda r: span is not None and r.sseqid == span[0] and span[1] <= min(r.sstart, r.send) and max(r.sstart, r.send) <= span[2], axis=1) if len(hx) else pd.Series(dtype=bool)
        hx = hx.assign(target=is_t.values if len(hx) else [])
        if g == "langdon": res["target_langdon"] = bool(hx[hx.role.isin(["X", "Y"])].target.any() and hx[hx.role == "C"].target.any())
        off = hx[~hx.target]
        res["off_hits"] += len(off)
        # посторонние продукты: аллельный + общий на одной последовательности навстречу, ≤ 1000 п.н.
        A_, C_ = off[off.role.isin(["X", "Y"])], off[off.role == "C"]
        for _, a in A_.iterrows():
            c = C_[(C_.sseqid == a.sseqid) & (C_.strand != a.strand)]
            if len(c) and any(abs(a.sstart - x) <= 1000 and ((a.strand == "+" and x > a.sstart) or (a.strand == "-" and x < a.sstart)) for x in c.sstart):
                res["off_products"] += 1
        for role in ["X", "Y", "C"]:
            hr = hx[hx.role == role]
            rows_spec.append(dict(маркер=d.marker, вариант_праймеров=d.design, праймер={"X": "FAM-аллель", "Y": "HEX-аллель", "C": "общий"}[role],
                                  геном=GEN[g], попаданий=len(hr), целевых=int(hr.target.sum()), нецелевых=int((~hr.target).sum()),
                                  координаты="; ".join(f"{r.sseqid}:{r.sstart}-{r.send}({'цель' if r.target else 'вне цели'}, несовп. {r.mm})" for _, r in hr.head(12).iterrows())
                                  + (" …" if len(hr) > 12 else "")))
    eval_.append(res)
E = pd.DataFrame(eval_)
DE = D.merge(E, on=["marker", "design"]).merge(T[["marker", "cand_id", "target_rank", "type"]], on="marker")
DE["grp"] = np.where(DE.cand_id.fillna("") == "", DE.marker, DE.cand_id)
DE["ok"] = (DE.off_products == 0) & DE.target_langdon
chosen = DE.sort_values(["grp", "ok", "off_products", "off_hits", "target_rank", "design"], ascending=[True, False, True, True, True, True]).groupby("grp").head(1)

pg = pd.read_csv(K + "parents_gt.tsv", sep="\t"); pg.columns = [re.sub(r"^#?\[\d+\]", "", c).replace(":GT", "") for c in pg.columns]
rename = dict(l.split() for l in open(os.path.join(OUT, "tables/02_chr_rename.txt"), encoding="utf-8") if l.strip()); pg["chr"] = pg.CHROM.map(rename)


def gt_call(gt, ref, alt, X, Y):
    if gt in ("./.", ".", "nan") or pd.isna(gt): return "нет вызова GBS"
    al = {"0": ref, "1": alt}; a = [al.get(x, "?") for x in re.split(r"[/|]", gt)]
    lab = lambda b: "FAM" if b == X else ("HEX" if b == Y else f"другой ({b})")
    return f"{'/'.join(a)} ({lab(a[0])}{'' if a[0] == a[1] else ' + ' + lab(a[1])})"


out = []
for _, c in chosen.iterrows():
    t = T[T.marker == c.marker].iloc[0]
    X, Y = t.allele_langdon, t.allele_group2
    al = dict(x.split(":", 1) for x in str(t.alleles_assemblies).split("; ") if ":" in x) if isinstance(t.alleles_assemblies, str) else {}
    exp = {"Langdon": f"{X}/{X} (FAM)"}
    if t.type.startswith("а"):
        sv = al.get("Svevo"); exp["Svevo"] = f"{sv}/{sv} ({'FAM' if sv == X else 'HEX' if sv == Y else 'другой'})" if sv else "нет продукта (нет копии гена в Svevo)"
    else:
        sv = TH["svevo"]; sv = sv[(sv.qseqid == c.marker)].sort_values("bitscore", ascending=False)
        b = None
        if len(sv):
            r = sv.iloc[0]; q = int(r.qstart) - 1
            for a_, s_ in zip(r.qseq, r.sseq):
                if a_ != "-":
                    if q == int(t.snp_index): b = s_.upper(); break
                    q += 1
        exp["Svevo"] = f"{b}/{b} ({'FAM' if b == X else 'HEX' if b == Y else 'другой'}; по сборке Svevo v1)" if b else "нет данных"
    x = pg[(pg.chr == t.chr) & (pg.POS == int(t.pos_langdon))]
    for tp, nm in PAR + [("TP_MO_216", "SVEVO (GBS коллекции)")]:
        exp[nm] = gt_call(str(x.iloc[0][tp]), x.iloc[0].REF, x.iloc[0].ALT, X, Y) if len(x) else "нет данных GBS (SNP маркера не в VCF)"
    out.append(dict(маркер=c.marker, тип=t.type, ген=t.gene if isinstance(t.gene, str) else "", вариант=(t.variant_note if isinstance(t.variant_note, str) and t.variant_note else
                    f"{t.rel} {X}>{Y}" if isinstance(t.rel, str) and t.rel else f"SNP {X}>{Y}"),
                    группы=(f"FAM ({X}): {t.group1}; HEX ({Y}): {t.group2}" if t.type.startswith("а") else ""),
                    зона=t.zone, координаты_Langdon=f"{t.chr}:{int(t.pos_langdon)}", аллель_FAM=X, аллель_HEX=Y,
                    праймер_FAM=FAM + c.as_x, праймер_HEX=HEX + c.as_y, общий_праймер=c.common,
                    Tm_FAM=c.tm_x, Tm_HEX=c.tm_y, Tm_общий=c.tm_c, ампликон_пн=int(c.amplicon), ориентация_аллельных=("прямая" if c.orient == "F" else "обратная"),
                    посторонних_продуктов=int(c.off_products), нецелевых_попаданий=int(c.off_hits), целевой_продукт_Langdon="да" if c.target_langdon else "нет",
                    статус=("годен по in silico проверке: целевой продукт в Langdon, посторонних продуктов нет"
                            + (f"; одиночных нецелевых посадок праймеров {int(c.off_hits)}" if c.off_hits else "")) if c.ok
                           else "есть нецелевые продукты или нет целевого — проверить",
                    аллели_сборок=t.alleles_assemblies if isinstance(t.alleles_assemblies, str) else "", **{f"ожидаемо_{k}": v for k, v in exp.items()}))
# гены, для которых KASP не подобран / не применим
for g, x in T[T.type.str.startswith("а")].groupby("cand_id"):
    if g in set(chosen.cand_id.dropna()): continue
    st = x.status.iloc[0] if (x.status != "мишень выбрана").all() else "праймеры с требованиями не подобраны ни для одной из мишеней"
    out.append(dict(маркер="", тип="а: ген-кандидат", ген=x.gene.iloc[0], вариант=str(x.group2_diff.iloc[0]), статус=st))
for _, x in T[T.type.str.startswith("б")].iterrows():      # SNP GWAS без подобранных праймеров
    if x.marker not in set(chosen.marker):
        dd = pd.read_csv(K + "designs.tsv", sep="	"); st = dd[dd.marker == x.marker].status
        out.append(dict(маркер=x.marker, тип=x.type, ген="", вариант=f"SNP {x.allele_langdon}>{x.allele_group2}", координаты_Langdon=f"{x.chr}:{int(x.pos_langdon)}",
                        статус=st.iloc[0] if len(st) else "праймеры не подобраны"))
# доминантная ПЦР по делеции CBF-B12/B14/B15
pr = pd.read_csv(os.path.join(OUT, "tables/11_pearce2013_primers.csv"), encoding="utf-8-sig").set_index("gene")
ins = pd.read_csv(os.path.join(OUT, "tables/16_cbf_primers_insilico.csv"), encoding="utf-8-sig") if os.path.exists(os.path.join(OUT, "tables/16_cbf_primers_insilico.csv")) else pd.DataFrame()
fp = pd.read_csv(os.path.join(OUT, "tables/14_f2_parents_genotypes.csv"), encoding="utf-8-sig")
fbc = [c for c in fp.columns if c.startswith("Fr-B2 (класс)")]       # GBS-гаплотип Fr-B2 родителей (3 SNP), без вывода о продукте ПЦР
fb2 = {r["parent"]: "GBS Fr-B2 (5B:521642333/522742783/532482598): " + "/".join(str(r[c]) for c in fbc) for _, r in fp.iterrows()}
for pair, role in [("CBF-B12-pair 2", "мишень делеции"), ("CBF-B14-pair 1", "мишень делеции"), ("CBF-B15-pair 1", "мишень делеции"), ("CBF-B9-pair 1", "положительный контроль")]:
    e = ins[ins.primer_pair == pair] if len(ins) else ins
    prod = "; ".join(f"{r['genome']}: {str(r['sizes_bp']) + ' п.н.' if r['product_in_FrB2_region'] == 'да' else 'нет продукта в районе Fr-B2'}"
                     for _, r in e.iterrows()) if len(e) else "см. 16_cbf_primers_insilico.csv"     # r['product'], не r.product (метод Series)
    out.append(dict(маркер=f"PCR_{pair.replace(' ', '_')}", тип="доминантная ПЦР (делеция CBF-B12/B14/B15, Fr-B2)", ген=pair.split("-pair")[0], вариант=role,
                    праймер_FAM="", праймер_HEX="", общий_праймер=f"F {pr.loc[pair, 'forward_5_3']}; R {pr.loc[pair, 'reverse_5_3']}",
                    ампликон_пн=int(pr.loc[pair, "expected_product_bp"]), статус="Pearce et al. 2013 (Table S1); ожидаемые продукты по сборкам — in silico ПЦР шага 4",
                    аллели_сборок=prod, ожидаемо_Langdon=next((x.split(": ", 1)[1] for x in prod.split("; ") if x.startswith("Langdon:")), ""),
                    ожидаемо_Svevo=next((x.split(": ", 1)[1] for x in prod.split("; ") if x.startswith("Svevo:")), ""),
                    **{f"ожидаемо_{nm}": fb2.get(tp, "нет в VCF") for tp, nm in PAR}))
P = pd.DataFrame(out)
P.to_csv(os.path.join(OUT, "tables/16_kasp_primers.csv"), index=False, encoding="utf-8-sig")
S = pd.DataFrame(rows_spec)
S = S.merge(chosen[["marker", "design"]].rename(columns={"marker": "маркер", "design": "вариант_праймеров"}), on=["маркер", "вариант_праймеров"])
S.to_csv(os.path.join(OUT, "tables/16_kasp_specificity.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250)
print(P[["маркер", "ген", "ампликон_пн", "посторонних_продуктов", "нецелевых_попаданий", "статус"]].to_string())
