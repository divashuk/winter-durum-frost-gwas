# -*- coding: utf-8 -*-
"""Раздел 16, шаги 5б и 6. Варианты по выравниваниям MAFFT (genes16/aln) и оценка их функциональности.
Опора — копия Langdon (референс VCF); позиции даются в координатах копии Langdon: g.−N — N п.н. до начала CDS (промотор 1 кб = g.−1000…−1),
g.N — от начала CDS (тело гена, для генов с интронами — включая интроны), p.N — белок Langdon, c.N — CDS Langdon.
Типы вариантов: SNP, инсерция/делеция (< 50 п.н.), протяжённая делеция/инсерция (≥ 50 п.н.), миссенс, синонимичный, стоп, сдвиг рамки,
белковый индел (в рамке), потеря гена (нет копии в синтенном окне), изменение числа копий (CNV), мотивы промотора.
Эффект в CDS: синонимичный/миссенс — по кодонам, сопоставленным через белковое выравнивание (только копии без сдвига рамки);
стоп и сдвиг рамки — по miniprot (StopCodon, Frameshift; 16_04k_integrate.py).
Домены: InterPro (Ensembl Plants BioMart) белка CS (для CBF — паралог CS) или Svevo, перенесённые на белок Langdon глобальным
выравниванием (BLOSUM62); для CBF дополнительно сигнатуры PKK/RPAGRxKFxETRHP и DSAWR (поиск по белку Langdon, ≤ 1 замены).
Консервативность позиции: выравнивание белка с ортологами ячменя, риса, Brachypodium (Compara): «консервативна» — у всех
ортологов (≥ 2) остаток как у Langdon; «частично» — у части; «нет» — ни у одного.
Grantham — матрица Grantham (1974).
Промотор (1 кб до CDS): мотивы CRT/DRE CCGAC (и GTCGG), ABRE ACGTG (и CACGT), ICE-box CANNTG; потеря/появление — мотив есть в одной из
двух последовательностей в той же колонке выравнивания.
КЛАСС (правило):
  сильный  — потеря гена; стоп; сдвиг рамки; белковый индел или миссенс в сигнатуре CBF (PKK…ETRHP, DSAWR);
  умеренный— миссенс в домене InterPro или в консервативной позиции при Grantham ≥ 60; белковый индел вне сигнатур; CNV;
             потеря/появление CRT/DRE или ABRE в промоторе 1 кб; протяжённый индел (≥ 50 п.н.) в промоторе или в теле гена с интронами
             (у CBF тело гена = CDS, его инделы учтены на уровне белка);
  слабый   — прочие миссенс; синонимичные; изменения ICE-box; прочие SNP/инделы промотора, интронов и фланков;
  «не оценивается (сборка)» — у сборки нет копии в синтенном окне, но есть на неразмещённых контигах (это не потеря гена).
Группы гаплотипов: по каждому гену — (2) носители одного варианта наибольшего класса (при равенстве — того, что есть у большего
числа сборок; варианты с тем же набором носителей перечисляются вместе), (1) сборки с белком как у Langdon, кроме группы 2.
Родители F2 — по GBS (VCF коллекции), если SNP GBS лежит в запросе гена и виден в выравнивании; иначе «по GBS не видно».
Выход: tables/16_variants_functional.csv, tables/16_haplotype_groups.csv."""
import os, re, sys
import numpy as np, pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
S16 = "/work/scratch/s16/"
ALN = os.path.join(OUT, "genes16", "aln")
CODON = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
                                                 "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
GR = """A R 112|A N 111|A D 126|A C 195|A Q 91|A E 107|A G 60|A H 86|A I 94|A L 96|A K 106|A M 84|A F 113|A P 27|A S 99|A T 58|A W 148|A Y 112|A V 64|R N 86|R D 96|R C 180|R Q 43|R E 54|R G 125|R H 29|R I 97|R L 102|R K 26|R M 91|R F 97|R P 103|R S 110|R T 71|R W 101|R Y 77|R V 96|N D 23|N C 139|N Q 46|N E 42|N G 80|N H 68|N I 149|N L 153|N K 94|N M 142|N F 158|N P 91|N S 46|N T 65|N W 174|N Y 143|N V 133|D C 154|D Q 61|D E 45|D G 94|D H 81|D I 168|D L 172|D K 101|D M 160|D F 177|D P 108|D S 65|D T 85|D W 181|D Y 160|D V 152|C Q 154|C E 170|C G 159|C H 174|C I 198|C L 198|C K 202|C M 196|C F 205|C P 169|C S 112|C T 149|C W 215|C Y 194|C V 192|Q E 29|Q G 87|Q H 24|Q I 109|Q L 113|Q K 53|Q M 101|Q F 116|Q P 76|Q S 68|Q T 42|Q W 130|Q Y 99|Q V 96|E G 98|E H 40|E I 134|E L 138|E K 56|E M 126|E F 140|E P 93|E S 80|E T 65|E W 152|E Y 122|E V 121|G H 98|G I 135|G L 138|G K 127|G M 127|G F 153|G P 42|G S 56|G T 59|G W 184|G Y 147|G V 109|H I 94|H L 99|H K 32|H M 87|H F 100|H P 77|H S 89|H T 47|H W 115|H Y 83|H V 84|I L 5|I K 102|I M 10|I F 21|I P 95|I S 142|I T 89|I W 61|I Y 33|I V 29|L K 107|L M 15|L F 22|L P 98|L S 145|L T 92|L W 61|L Y 36|L V 32|K M 95|K F 102|K P 103|K S 121|K T 78|K W 110|K Y 85|K V 97|M F 28|M P 87|M S 135|M T 81|M W 67|M Y 36|M V 21|F P 114|F S 155|F T 103|F W 40|F Y 22|F V 50|P S 74|P T 38|P W 147|P Y 110|P V 68|S T 58|S W 177|S Y 144|S V 124|T W 128|T Y 92|T V 69|W Y 37|W V 88|Y V 55"""
GRANTHAM = {}
for x in GR.split("|"):
    a, b, v = x.split(); GRANTHAM[(a, b)] = GRANTHAM[(b, a)] = int(v)
MOTIFS = {"CRT/DRE": [r"CCGAC", r"GTCGG"], "ABRE": [r"ACGTG", r"CACGT"], "ICE-box": [r"CA..TG"]}
LABEL = {"ASM5281831v1": "PI 192051", "Zavitan_WEW2.1": "Zavitan", "Dicoccum_lpTriDico3": "dicoccum lpTriDico3", "CS_v1.0": "Chinese Spring",
         "Stanley": "CDC Stanley", "Landmark": "CDC Landmark", "SY_Mattis": "SY Mattis", "Norin61": "Norin 61"}   # как в 16_04k_integrate.py
SIG = {"сигнатура PKKxPAGRxKFxETRHP": "PKKRPAGRTKFKETRHP", "сигнатура DSAWR": "DSAWR"}


def read_aln(path):
    if not os.path.exists(path): return [], []
    names, seqs, cur = [], [], []
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            if names: seqs.append("".join(cur).upper())
            names.append(line[1:]); cur = []
        else: cur.append(line)
    if names: seqs.append("".join(cur).upper())
    return names, seqs


def ref_index(names):
    for i, n in enumerate(names):
        if n.startswith("Langdon|") and "вне_синтении" not in n: return i
    return None


def colmap(s):
    """колонка -> позиция в негэпованной последовательности (1-based) или None"""
    out, k = [], 0
    for ch in s:
        if ch != "-": k += 1; out.append(k)
        else: out.append(None)
    return out


def find_sig(prot, sig):
    """позиция сигнатуры в белке: ≤ 3 замены для PKK…ETRHP (17 а.о.), ≤ 1 для DSAWR"""
    maxmm = 3 if len(sig) > 10 else 1; best = None
    for i in range(len(prot) - len(sig) + 1):
        mm = sum(1 for a, b in zip(prot[i:i + len(sig)], sig) if a != b)
        if mm <= maxmm and (best is None or mm < best[1]): best = (i + 1, mm)
    return (best[0], best[0] + len(sig) - 1) if best else None


def pair_map(a, b):
    """глобальное выравнивание белков a -> b: словарь позиция_a -> позиция_b (1-based)"""
    al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -0.5; al.mode = "global"
    aln = al.align(a, b)[0]; m = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for k in range(e1 - s1): m[s1 + k + 1] = s2 + k + 1
    return m


def short(n):     # метка сборки как в 16_gene_copies_by_genome.csv (пробелы вместо «_» из заголовков fasta)
    p = n.split("|"); return p[0].replace("_", " ") + ("" if p[2] == "копия1" else "/" + p[2].replace("копия", "к")) + ("*вне" if "вне_синтении" in n else "")


def genome_of(n): return n.split("|")[0].replace("_", " ")


def main():
    names = pd.read_csv(S16 + "q/candidate_names.tsv", sep="\t")
    fb = pd.read_csv(S16 + "q/file_base.tsv", sep="\t", header=None, names=["base", "cand_id"])
    A = pd.read_csv(S16 + "copies/all_copies.tsv", sep="\t")
    T = pd.read_csv(os.path.join(OUT, "tables/16_gene_copies_by_genome.csv"))
    cand = pd.read_csv(S16 + "q/candidates.tsv", sep="\t").set_index("cand_id")
    ipr_cs = pd.read_csv(S16 + "biomart/cs_interpro.tsv", sep="\t").drop_duplicates(["ensembl_gene_id", "interpro", "interpro_start"])
    ipr_sv = pd.read_csv(S16 + "biomart/svevo_interpro.tsv", sep="\t").drop_duplicates(["ensembl_gene_id", "interpro", "interpro_start"])
    cs_pep = {re.search(r"gene:(\S+)", h).group(1): s for h, s in read_fa_full(S16 + "seq/cs_pep.fa").items()}
    sv_pep = {h.split()[0]: s for h, s in read_fa_full(S16 + "seq/svevo_pep.fa").items()}
    pg = pd.read_csv(S16 + "kasp/parents_gt.tsv", sep="\t")
    pg.columns = [re.sub(r"^#?\[\d+\]", "", c).replace(":GT", "") for c in pg.columns]
    PNAME = {"TP_MO_081": "Цель", "TP_MO_044": "Андромеда", "TP_MO_054": "Прибуткова", "TP_MO_014": "Синьора", "TP_MO_005": "Алена",
             "TP_MO_037": "Степной янтарь", "TP_MO_171": "ANTALIS", "TP_MO_190": "PAPRICA", "TP_MO_196": "CORTEZ", "TP_MO_302": "PALEOTTO",
             "TP_MO_304": "L2870 VILLOSA", "TP_MO_350": "RAMIRES", "TP_MO_216": "SVEVO (коллекция)"}
    rename = dict(l.split() for l in open(os.path.join(OUT, "tables/02_chr_rename.txt"), encoding="utf-8") if l.strip())
    pg["chr"] = pg.CHROM.map(rename)
    V, G = [], []
    for _, fbr in fb.iterrows():
        base, cid = fbr.base, fbr.cand_id
        nm = names[names.cand_id == cid].iloc[0]; gene = nm.ген; c = cand.loc[cid]
        is_cbf = cid.startswith("CBF_")
        pn, ps = read_aln(os.path.join(ALN, f"{base}_protein.aln.fa"))
        gn, gs = read_aln(os.path.join(ALN, f"{base}_genomic.aln.fa"))
        on, os_ = read_aln(os.path.join(ALN, f"{base}_protein_orth.aln.fa"))
        ri = ref_index(pn); gri = ref_index(gn)
        if ri is None or gri is None:
            print("нет копии Langdon:", gene); continue
        rprot = ps[ri].replace("-", "")
        # ---------- домены на белке Langdon ----------
        dom = []
        src_id = nm.ID_CS_RefSeq_v1_1 if is_cbf else nm.ID_Svevo
        src_pep = cs_pep.get(src_id) if is_cbf else sv_pep.get(src_id)
        src_ipr = ipr_cs[ipr_cs.ensembl_gene_id == src_id] if is_cbf else ipr_sv[ipr_sv.ensembl_gene_id == src_id]
        if src_pep:
            m = pair_map(src_pep.rstrip("*"), rprot)
            for _, d in src_ipr.dropna(subset=["interpro_start"]).iterrows():
                a, b = m.get(int(d.interpro_start)), m.get(int(d.interpro_end))
                if a and b and (b - a + 1) < 0.8 * len(rprot):     # записи на весь белок (семейство, суперсемейство) — не домен
                    dom.append((f"{d.interpro_description} ({d.interpro})", a, b))
        if is_cbf:
            for k, sig in SIG.items():
                x = find_sig(rprot, sig)
                if x: dom.append((k, x[0], x[1]))
        def dom_at(p):
            return "; ".join(sorted({d[0] for d in dom if d[1] <= p <= d[2]}))
        # ---------- консервативность ----------
        cons = {}
        if on:
            oi = ref_index(on)
            if oi is not None:
                cm = colmap(os_[oi]); orth = [j for j, n in enumerate(on) if "|" in n and n.split("|")[0] in ("hordeum_vulgare", "oryza_sativa", "brachypodium_distachyon")]
                for col, p in enumerate(cm):
                    if p is None: continue
                    res = [os_[j][col] for j in orth if os_[j][col] != "-"]
                    if len(res) >= 2:
                        same = sum(r == os_[oi][col] for r in res)
                        cons[p] = "консервативна" if same == len(res) else ("частично" if same else "нет")
        # ---------- белок: миссенс, инделы, стоп ----------
        rcm = colmap(ps[ri]); events = {}
        for j, (n, s) in enumerate(zip(pn, ps)):
            if j == ri: continue
            first = next((k for k, ch in enumerate(s) if ch != "-"), None); last = len(s) - next((k for k, ch in enumerate(s[::-1]) if ch != "-"), 0)
            k = first or 0
            while k < last:
                r, a = ps[ri][k], s[k]
                if r != "-" and a != "-" and r != a:
                    key = ("missense" if a != "*" else "stop", rcm[k], r, a); events.setdefault(key, set()).add(n); k += 1
                elif r == "-" and a != "-" or r != "-" and a == "-":
                    t = "ins" if r == "-" else "del"; k0 = k
                    while k < last and ((ps[ri][k] == "-") == (t == "ins")) and ((s[k] == "-") == (t == "del")): k += 1
                    pos = next((rcm[x] for x in range(k0, len(rcm)) if rcm[x]), len(rprot))
                    key = (f"белковый {('индел: вставка' if t == 'ins' else 'индел: делеция')}", pos, ps[ri][k0:k].replace("-", ""), s[k0:k].replace("-", ""))
                    events.setdefault(key, set()).add(n)
                else: k += 1
        for (t, pos, r, a), carriers in events.items():
            if t in ("missense", "stop"):
                eff = "миссенс" if t == "missense" else "стоп"
                gr = GRANTHAM.get((r, a)) if t == "missense" else None
                d = dom_at(pos); cv = cons.get(pos, "нет данных")
                if t == "stop" or "сигнатура" in d: cls = "сильный"
                elif d or (cv == "консервативна" and (gr or 0) >= 60): cls = "умеренный"
                else: cls = "слабый"
                V.append(dict(ген=gene, cand_id=cid, уровень="белок", вариант=f"p.{r}{pos}{'*' if a == '*' else a}", позиция_Langdon=f"p.{pos}",
                              эффект=eff, домен=d, консервативность=cv, Grantham=gr, мотив="", класс=cls,
                              сборки=", ".join(sorted(short(x) for x in carriers)), n_сборок=len({genome_of(x) for x in carriers})))
            else:
                d = "; ".join(sorted({dom_at(p) for p in range(pos, pos + max(1, len(r)))} - {""}))
                cls = "сильный" if "сигнатура" in d else "умеренный"
                V.append(dict(ген=gene, cand_id=cid, уровень="белок", вариант=f"{t} {len(r) or len(a)} а.о. ({r or a})", позиция_Langdon=f"p.{pos}",
                              эффект=t, домен=d, консервативность=cons.get(pos, "нет данных"), Grantham=None, мотив="", класс=cls,
                              сборки=", ".join(sorted(short(x) for x in carriers)), n_сборок=len({genome_of(x) for x in carriers})))
        # ---------- синонимичные (кодоны через белковое выравнивание) ----------
        cdsd = {}
        for _, r in A[A.cand_id == cid].iterrows():
            lab = r.genome
            cdsd.setdefault(lab, []).append(r)
        cn, cs_ = read_aln(os.path.join(ALN, f"{base}_cds.aln.fa"))
        cdseq = {n: s.replace("-", "") for n, s in zip(cn, cs_)}
        rc_ = cdseq.get(pn[ri], "")
        syn = {}
        ref_fs = A[(A.cand_id == cid) & (A.genome == "Langdon")].orf_why.fillna("").str.contains("сдвиг рамки").any()
        if len(rc_) >= 3 * len(rprot) and len(rc_) <= 3 * len(rprot) + 3 and not ref_fs:   # у Langdon сдвиг рамки — синонимичные не оцениваются
            for j, (n, s) in enumerate(zip(pn, ps)):
                if j == ri: continue
                ac = cdseq.get(n, ""); ap = s.replace("-", "")
                if len(ac) < 3 * len(ap) or len(ac) > 3 * len(ap) + 3: continue     # сдвиг рамки / несоответствие — пропуск
                acm = colmap(s)
                for k in range(len(s)):
                    if rcm[k] and acm[k] and ps[ri][k] == s[k]:
                        x, y = rc_[3 * (rcm[k] - 1):3 * rcm[k]], ac[3 * (acm[k] - 1):3 * acm[k]]
                        if x != y and len(x) == len(y) == 3:
                            for q in range(3):
                                if x[q] != y[q]: syn.setdefault((3 * (rcm[k] - 1) + q + 1, x[q], y[q], ps[ri][k]), set()).add(n)
        for (cpos, rb, ab, aa), carriers in syn.items():
            V.append(dict(ген=gene, cand_id=cid, уровень="CDS", вариант=f"c.{cpos}{rb}>{ab} (p.{aa}{(cpos - 1)//3 + 1}=)", позиция_Langdon=f"c.{cpos}",
                          эффект="синонимичный", домен=dom_at((cpos - 1) // 3 + 1), консервативность="", Grantham=None, мотив="", класс="слабый",
                          сборки=", ".join(sorted(short(x) for x in carriers)), n_сборок=len({genome_of(x) for x in carriers})))
        # ---------- геномные: SNP/инделы вне CDS, промоторные мотивы ----------
        L = int(c.end - c.start + 1)
        gcm = colmap(gs[gri]); gev = {}
        def rel(i):          # позиция в копии Langdon -> обозначение
            return f"g.{i - 2001}" if i <= 2000 else (f"g.{i - 2000}" if i <= 2000 + L else f"g.*{i - 2000 - L}")
        def where(i):
            return "промотор (1 кб)" if 1000 < i <= 2000 else ("5′-фланк (1–2 кб)" if i <= 1000 else ("тело гена" if i <= 2000 + L else "3′-фланк"))
        for j, (n, s) in enumerate(zip(gn, gs)):
            if j == gri: continue
            first = next((k for k, ch in enumerate(s) if ch != "-"), 0); last = len(s) - next((k for k, ch in enumerate(s[::-1]) if ch != "-"), 0)
            k = first
            while k < last:
                r, a = gs[gri][k], s[k]
                if r != "-" and a != "-":
                    if r != a and "N" not in (r, a): gev.setdefault(("SNP", gcm[k], r, a), set()).add(n)
                    k += 1
                elif r == "-" and a == "-": k += 1
                else:
                    t = "ins" if r == "-" else "del"; k0 = k
                    while k < last and ((gs[gri][k] == "-") == (t == "ins")) and ((s[k] == "-") == (t == "del")) and not (gs[gri][k] == "-" and s[k] == "-"): k += 1
                    pos = next((gcm[x] for x in range(k0, len(gcm)) if gcm[x]), None) or len(gs[gri].replace("-", ""))
                    gev.setdefault((t, pos, k - k0), set()).add(n)
        for key, carriers in gev.items():
            t, pos = key[0], key[1]
            w = where(pos)
            if w == "тело гена" and is_cbf: continue          # у CBF (без интронов) тело гена = CDS — эффекты учтены на уровне белка/CDS
            if t == "SNP":
                var, eff = f"{rel(pos)}{key[2]}>{key[3]}", "SNP"
            else:
                ln = key[2]; eff = ("протяжённая " if ln >= 50 else "") + ("инсерция" if t == "ins" else "делеция")
                var = f"{rel(pos)} {eff} {ln} п.н."
            if w == "тело гена" and not is_cbf and t == "SNP": continue      # тело многоэкзонного гена: эффект в CDS оценён по белку
            if "протяжённая" in eff and w in ("промотор (1 кб)", "тело гена"):
                cls = "умеренный"        # у генов с интронами тело гена включает интроны: эффект на белок — по белковому выравниванию
            else:
                cls = "слабый"
            V.append(dict(ген=gene, cand_id=cid, уровень="геномная", вариант=var, позиция_Langdon=rel(pos), эффект=eff + f" ({w})", домен="",
                          консервативность="", Grantham=None, мотив="", класс=cls, сборки=", ".join(sorted(short(x) for x in carriers)),
                          n_сборок=len({genome_of(x) for x in carriers})))
        # мотивы промотора
        def motifs(s):
            u = s.replace("-", ""); cm = [k for k, ch in enumerate(s) if ch != "-"]; out = set()
            for mname, pats in MOTIFS.items():
                for p in pats:
                    for m_ in re.finditer(f"(?=({p}))", u):
                        i = m_.start() + 1
                        if 1000 < i <= 2000: out.add((mname, cm[i - 1]))
            return out
        rm = motifs(gs[gri]); mev = {}
        for j, (n, s) in enumerate(zip(gn, gs)):
            if j == gri: continue
            first = next((k for k, ch in enumerate(s) if ch != "-"), 0)
            am = motifs(s) if colmap(s)[first] is not None else set()
            for x in rm - am:
                if s[x[1]] != "-" or True: mev.setdefault(("потеря", x), set()).add(n)
            for x in am - rm: mev.setdefault(("появление", x), set()).add(n)
        for (t, (mname, col)), carriers in mev.items():
            # проверка покрытия: у сборки колонка должна быть внутри её выровненной части
            carriers = {n for n in carriers if any(ch != "-" for ch in gs[gn.index(n)][:col + 1]) and any(ch != "-" for ch in gs[gn.index(n)][col:])}
            if not carriers: continue
            p = next((gcm[x] for x in range(col, len(gcm)) if gcm[x]), None)
            if p is None or not (1000 < p <= 2000): continue     # только промотор 1 кб в координатах Langdon
            V.append(dict(ген=gene, cand_id=cid, уровень="промотор", вариант=f"{t} {mname} ({rel(p) if p else '—'})", позиция_Langdon=rel(p) if p else "",
                          эффект=f"{t} мотива {mname}", домен="", консервативность="", Grantham=None, мотив=mname,
                          класс="умеренный" if mname in ("CRT/DRE", "ABRE") else "слабый", сборки=", ".join(sorted(short(x) for x in carriers)),
                          n_сборок=len({genome_of(x) for x in carriers})))
        # ---------- потеря гена, CNV, сдвиг рамки ----------
        tt = T[T.cand_id == cid]
        nL = int(tt[tt.сборка == "Langdon"].копий_в_синтенном_окне.iloc[0]) if (tt.сборка == "Langdon").any() else 1
        unpl_only = tt[(tt.копий_в_синтенном_окне == 0) & (tt.копий_на_неразмещённых_контигах > 0)].сборка.tolist()
        lost = tt[(tt.копий_в_синтенном_окне == 0) & (tt.копий_на_неразмещённых_контигах == 0)].сборка.tolist()
        if unpl_only:      # копия есть, но не в псевдомолекуле: положение неизвестно — не «потеря гена», класс не присваивается
            V.append(dict(ген=gene, cand_id=cid, уровень="копии", вариант="копии только на неразмещённых контигах", позиция_Langdon="",
                          эффект="положение в сборке не определено", домен="", консервативность="", Grantham=None, мотив="",
                          класс="не оценивается (сборка)", сборки=", ".join(unpl_only), n_сборок=len(unpl_only)))
        cnv = tt[(tt.копий_в_синтенном_окне > 0) & (tt.копий_в_синтенном_окне != nL)]
        if lost:
            V.append(dict(ген=gene, cand_id=cid, уровень="копии", вариант="нет копии в синтенном окне", позиция_Langdon="", эффект="потеря гена", домен="",
                          консервативность="", Grantham=None, мотив="", класс="сильный", сборки=", ".join(lost), n_сборок=len(lost)))
        for _, r in cnv.iterrows():
            V.append(dict(ген=gene, cand_id=cid, уровень="копии", вариант=f"{r.копий_в_синтенном_окне} копий (Langdon {nL})", позиция_Langdon="",
                          эффект="изменение числа копий", домен="", консервативность="", Grantham=None, мотив="", класс="умеренный", сборки=r.сборка, n_сборок=1))
        fsx = A[(A.cand_id == cid) & A.orf_why.fillna("").str.contains("сдвиг рамки") & (A.dna_ident >= 0.95)]
        if ref_fs: fsx = fsx.iloc[0:0]          # сдвиг рамки есть и у Langdon — это не различие между сборками
        if len(fsx):
            V.append(dict(ген=gene, cand_id=cid, уровень="CDS", вариант="сдвиг рамки (miniprot)", позиция_Langdon="", эффект="сдвиг рамки", домен="",
                          консервативность="", Grantham=None, мотив="", класс="сильный",
                          сборки=", ".join(sorted({LABEL.get(x, x) for x in fsx.genome})), n_сборок=fsx.genome.nunique()))
        # ---------- группы гаплотипов ----------
        hap = {}
        for n, s in zip(pn, ps): hap.setdefault(s.replace("-", ""), []).append(n)
        refhap = hap[rprot]
        v = pd.DataFrame([x for x in V if x["cand_id"] == cid])
        rank = {"сильный": 3, "умеренный": 2, "слабый": 1}
        if len(v):
            v["r"] = v.класс.map(rank); top = v.r.max()
            tv = v[v.r == top].sort_values("n_сборок", ascending=False, kind="stable")
            first = tv.iloc[0]          # группа 2 — носители одного сильнейшего варианта (при равенстве — у большего числа сборок)
            g2 = sorted({x.split("/")[0] for x in str(first.сборки).replace("*вне", "").split(", ") if x})
            same = tv[tv.сборки == first.сборки]
            what = "; ".join(same.вариант.head(6)) + (" …" if len(same) > 6 else "") + (f" (ещё вариантов этого класса у других сборок: {len(tv) - len(same)})" if len(tv) > len(same) else "")
        else:
            top, g2, what = 0, [], "различий нет"
        g1 = sorted({genome_of(x) for x in refhap} - set(g2))
        # родители F2 по GBS
        chrn = c.chr; qs, qe = int(c.qstart), int(c.qend)
        x = pg[(pg.chr == chrn) & pg.POS.between(qs, qe)]
        par = []
        for _, s in x.iterrows():
            al = {PNAME[k]: s[k] for k in PNAME if k in s.index}
            par.append(f"{chrn}:{s.POS} {s.REF}/{s.ALT}: " + ", ".join(f"{k} {v_}" for k, v_ in al.items()))
        G.append(dict(ген=gene, cand_id=cid, приоритет=nm.приоритет, макс_класс={3: "сильный", 2: "умеренный", 1: "слабый", 0: "—"}[top],
                      группа_1_как_Langdon=", ".join(g1), группа_2=", ".join(g2), отличие_группы_2=what,
                      гаплотипов_белка=len(hap), родители_F2_GBS="; ".join(par) if par else "по GBS не видно (SNP GBS в гене нет)"))
        print(gene, "вариантов", len(v), "гаплотипов белка", len(hap))
    Vd = pd.DataFrame(V); Gd = pd.DataFrame(G)
    Vd.to_csv(os.path.join(OUT, "tables/16_variants_functional.csv"), index=False, encoding="utf-8-sig")
    Gd.to_csv(os.path.join(OUT, "tables/16_haplotype_groups.csv"), index=False, encoding="utf-8-sig")
    return Vd, Gd


def read_fa_full(path):
    d, n = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"): n = line[1:].strip(); d[n] = []
        elif n: d[n].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}


if __name__ == "__main__":
    main()
