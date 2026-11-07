# -*- coding: utf-8 -*-
"""Раздел 17, шаг 7в. Мишени KASP и шаблоны (Langdon ± 150 п.н.) -> s17/kasp/targets.tsv, templates.fa.
(а) ген-кандидат — scripts/16_07b_kasp_targets.py без изменения метода (колонка геномного выравнивания, полностью разделяющая группы 1 и 2
    tables/17_haplotype_groups.csv; маска — переменные колонки сборок) для генов части A; TRITD1Bv1G205280 (основной ген 1B; в разделе 16 имел
    приоритет 3 и KASP не подбирался) — по выравниваниям и группам раздела 16.
(б) SNP GWAS — в 16_07b список был задан для 5B; здесь — s17/kasp/snps_b.tsv (17_11a_snp_info.R): ведущие SNP локусов части A (7B — 7B:348683836
    и 7B:233439923) и 5B:535618213 (повтор с Tm 57–63 °C и ампликоном ≤ 150 п.н. — 17_14c). Заместитель 5B:535618213 с r² ≥ 0.8 в IT —
    только 5B:535723924; KASP на него подобран в разделе 16 (KASP_5B_535723924) и повторно не подбирается.
    Маска: все SNP обоих исходных VCF в ± 150 п.н. (s17/kasp/gwas_flank_vcf_snps.tsv; в 16 таких SNP не было, маска была пустой).
(в) покрытие гена (для части D: «на каждый ген-кандидат — не меньше одного маркера») — гены разделов 16 и 17 (приоритет 1–2, без фрагментов
    и повтора 513 Мб), у которых нет мишени типа (а) и KASP раздела 16: колонка геномного выравнивания (только синтенные копии), где у
    Langdon нет гэпа и среди сборок ровно два нуклеотида; порядок предпочтения — полиморфизм среди 6 сборок твёрдой пшеницы, тело гена /
    CDS, затем промотор 1 кб, затем фланки; без гэпов в ± 10 колонках; больше носителей минорного нуклеотида; ближе к середине гена;
    до 5 мишеней на ген не ближе 20 п.н. Маска — как в (а). Это не обязательно функциональный вариант.
Праймеры для всех типов подбираются одним кодом 16_07d (17_14c)."""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
S16 = lib.S16
NEWB = '''# (б) SNP GWAS — список раздела 17
SB = pd.read_csv(S16 + "kasp/snps_b.tsv", sep="\\t")
SB = SB[~SB.role.str.startswith("заместитель")]
pg = pd.read_csv(S16 + "kasp/parents_gt.tsv", sep="\\t")
pg.columns = [re.sub(r"^#?\\[\\d+\\]", "", c).replace(":GT", "") for c in pg.columns]
rename = dict(l.split() for l in open(os.path.join(OUT, "tables/02_chr_rename.txt"), encoding="utf-8") if l.strip())
pg["chr"] = pg.CHROM.map(rename)
vs = pd.read_csv(S16 + "kasp/gwas_flank_vcf_snps.tsv", sep="\\t", header=None, names=["CHROM", "POS", "REF", "ALT", "vcf"]); vs["chr"] = vs.CHROM.map(rename)
fl = {}
name = None
for line in open(S16 + "kasp/gwas_flanks.fa", encoding="utf-8"):
    if line.startswith(">"): name = line[1:].strip(); fl[name] = ""
    else: fl[name] += line.strip().upper()
for _, s in SB.iterrows():
    r = pg[(pg.chr == s.chr) & (pg.POS == int(s.pos))].iloc[0]
    t = fl[f"SNP_{s.chr}_{int(s.pos)}"]
    assert t[FL] == r.REF, (s.snp, t[FL], r.REF)
    lo = int(s.pos) - FL
    masked = sorted({int(p) - lo for p in vs[(vs.chr == s.chr) & vs.POS.between(lo, int(s.pos) + FL) & (vs.POS != int(s.pos))].POS})
    mk = f"KASP_{s.chr}_{int(s.pos)}"
    rows.append(dict(marker=mk, type="б: SNP GWAS", gene="", cand_id="", target_rank=1, status="мишень выбрана", chr=s.chr, pos_langdon=int(s.pos), rel="",
                     zone="SNP GBS", allele_langdon=r.REF, allele_group2=r.ALT, variant_note=f"SNP GBS коллекции (FINAL-181124); {s.role}; локус {s.locus}",
                     template_strand="+", snp_index=FL, masked=",".join(map(str, masked))))
    tmpl.append((mk, t))
# (в) покрытие гена
K16 = pd.read_csv(os.path.join(OUT, "tables/16_kasp_primers.csv"))
have16 = set(K16[K16.маркер.notna() & K16.тип.astype(str).str.startswith("а")].ген)
have_a = {r["cand_id"] for r in rows if r.get("status") == "мишень выбрана" and str(r.get("type", "")).startswith("а")}
DUR = {"Langdon", "Svevo", "Kronos", "PI 192051", "Chili", "Mahmoudi"}
need = []
for sec, SS, AL_ in [("16", "@S16@", os.path.join(OUT, "genes16", "aln")), ("17", S16, v.ALN)]:
    nn = pd.read_csv(SS + "q/candidate_names.tsv", sep="\\t"); ff = pd.read_csv(SS + "q/file_base.tsv", sep="\\t", header=None, names=["base", "cand_id"])
    cc = pd.read_csv(SS + "q/candidates.tsv", sep="\\t").set_index("cand_id")
    for _, f in ff.iterrows():
        nm = nn[nn.cand_id == f.cand_id]
        if nm.empty: continue
        nm = nm.iloc[0]
        if not str(nm.приоритет).startswith(("1", "2")) or "повтор" in str(nm.приоритет): continue
        if sec == "16" and (nm.ген in have16 or f.cand_id == "TRITD1Bv1G205280"): continue
        if sec == "17" and f.cand_id in have_a: continue
        need.append((sec, f.base, f.cand_id, nm, cc.loc[f.cand_id], AL_))
for sec, base, cid, nm, c, AL_ in need:
    gn, gs = v.read_aln(os.path.join(AL_, f"{base}_genomic.aln.fa"))
    ri = v.ref_index(gn) if gn else None
    rec = dict(marker="", type="в: покрытие гена (SNP между сборками)", gene=nm.ген, cand_id=cid, group1="", group2="", group2_diff=f"раздел {sec}")
    if ri is None or len(gn) < 2:
        rows.append({**rec, "status": "нет выравнивания с другими сборками — мишени нет"}); continue
    lab = [n.split("|")[0].replace("_", " ") for n in gn]
    idur = [j for j, l in enumerate(lab) if l in DUR]
    L = int(c.end - c.start + 1); refc = v.colmap(gs[ri]); cands = []
    var_cols = [k for k in range(len(gs[ri])) if len({s[k] for s in gs if s[k] != "-"}) > 1 or any(s[k] == "-" for s in gs)]
    for k in range(len(gs[ri])):
        col = [s[k] for s in gs]
        if col[ri] in "-N": continue
        al = [x for x in col if x not in "-N"]
        if len(set(al)) != 2: continue
        x = col[ri]; y = next(b for b in set(al) if b != x)
        near_gap = any(any(s[q] == "-" for s in gs) for q in range(max(0, k - 10), min(len(gs[ri]), k + 11)))
        dpoly = len({gs[j][k] for j in idur if gs[j][k] not in "-N"}) == 2
        p = refc[k]; zone = 0 if 2000 < p <= 2000 + L else (1 if 1000 < p <= 2000 else 2)
        cands.append(((0 if dpoly else 1, zone, near_gap, -min(al.count(x), al.count(y)), abs(p - (2000 + L / 2))), k, p, x, y, dpoly))
    if not cands:
        rows.append({**rec, "status": "нет двуаллельной колонки без гэпа у Langdon"}); continue
    picked = []
    for cnd in sorted(cands):
        if all(abs(cnd[2] - q[2]) >= 20 for q in picked): picked.append(cnd)
        if len(picked) == 5: break
    ref = gs[ri].replace("-", "")
    for rank, (_, k, p, x, y, dpoly) in enumerate(picked, 1):
        lo, hi = max(0, p - 1 - FL), min(len(ref), p + FL)
        t = ref[lo:hi]; snp_i = p - 1 - lo
        masked = sorted({refc[q] - 1 - lo for q in var_cols if refc[q] and q != k and lo <= refc[q] - 1 < hi})
        rel = f"g.{p - 2001}" if p <= 2000 else (f"g.{p - 2000}" if p <= 2000 + L else f"g.*{p - 2000 - L}")
        lc = int(c.qstart) + p - 1 if c.strand == "+" else int(c.qend) - p + 1
        zn = "CDS / тело гена" if 2000 < p <= 2000 + L else ("промотор 1 кб" if 1000 < p <= 2000 else "фланк")
        g1 = sorted({lab[j] for j in range(len(gs)) if gs[j][k] == x}); g2 = sorted({lab[j] for j in range(len(gs)) if gs[j][k] == y})
        mk = f"KASP_{re.sub(r'[^A-Za-z0-9]+', '_', nm.ген)[:28]}_{c.chr}_{lc}"
        rows.append({**rec, "marker": mk, "target_rank": rank, "status": "мишень выбрана", "chr": c.chr, "pos_langdon": lc, "rel": rel, "zone": zn,
                     "group1": ", ".join(g1), "group2": ", ".join(g2), "allele_langdon": x, "allele_group2": y,
                     "variant_note": f"{rel}{x}>{y}; полиморфизм среди твёрдой пшеницы: {'да' if dpoly else 'нет'}",
                     "alleles_assemblies": "; ".join(f"{a}:{gs[j][k]}" for j, a in enumerate(lab)), "template_strand": "ген", "snp_index": snp_i,
                     "masked": ",".join(map(str, masked))})
        tmpl.append((mk, t))
'''.replace("@S16@", S16)
src = lib.source("16_07b_kasp_targets.py")
old_b = src[src.index("# (б) SNP GWAS"):src.index("T = pd.DataFrame(rows)")]
add205 = (lambda what, extra: f'pd.concat([{what}, {extra}])')
lib.run16("16_07b_kasp_targets.py", [
    ('S16 = "/work/scratch/s16/"', f'S16 = "{lib.S17}"'),
    (lib.LOAD_V, "pass"),
    ('names = pd.read_csv(S16 + "q/candidate_names.tsv", sep="\\t")',
     f'names = pd.concat([pd.read_csv(S16 + "q/candidate_names.tsv", sep="\\t"), pd.read_csv("{S16}q/candidate_names.tsv", sep="\\t")'
     '.query("cand_id == \'TRITD1Bv1G205280\'").assign(приоритет="1 (основной ген 1B в разделе 17)")])'),
    ('fb = pd.read_csv(S16 + "q/file_base.tsv", sep="\\t", header=None, names=["base", "cand_id"])',
     'fb = pd.concat([pd.read_csv(S16 + "q/file_base.tsv", sep="\\t", header=None, names=["base", "cand_id"]), '
     'pd.DataFrame([["TRITD1Bv1G205280", "TRITD1Bv1G205280"]], columns=["base", "cand_id"])])'),
    ('cand = pd.read_csv(S16 + "q/candidates.tsv", sep="\\t").set_index("cand_id")',
     f'cand = pd.concat([pd.read_csv(S16 + "q/candidates.tsv", sep="\\t"), pd.read_csv("{S16}q/candidates.tsv", sep="\\t")'
     '.query("cand_id == \'TRITD1Bv1G205280\'")]).set_index("cand_id")'),
    ('os.path.join(OUT, "tables/16_haplotype_groups.csv")', 'os.path.join(OUT, "tables/17_haplotype_groups.csv")'),
    ('os.path.join(OUT, "tables/16_variants_functional.csv")', 'os.path.join(OUT, "tables/17_variants_functional.csv")'),
    (old_b, NEWB)], inject={"v": lib.variants17()})
