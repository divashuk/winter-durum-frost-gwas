# -*- coding: utf-8 -*-
"""Раздел 16, шаг 4л. Сводка копий генов-кандидатов по сборкам и последовательности genes16/.
Вход по каждой сборке g (сервер -> WSL s16/p2, s16/copies):
  copies/<g>_hits.tsv            — копии по ДНК (16_04f_copies.py: синтения, идентичность, покрытие тела гена);
  p2/<g>_win.gff, _win_long.gff  — miniprot (--gff --trans) белков CBF (-G 3k) и прочих кандидатов (-G 20k) на синтенных окнах;
  p2/<g>_win.fa, _copies.fa      — окна и последовательности копий (ген + 2 кб + 1 кб в ориентации гена).
Правила:
  * белок и CDS копии — лучшее (по счёту miniprot) выравнивание, перекрывающее тело гена копии ≥ 50 %; для CBF — любой белок CBF,
    для прочих — белок Svevo того же гена, иначе ортолог CS; «ORF полная» — покрытие белка ≥ 0.90, без сдвигов рамки и стопов;
  * число копий паралога в сборке — копии в синтенном окне (16_04f); «вне синтении» — отдельно;
  * у мягкой пшеницы гомеолог D: попадание запроса на хромосому D с идентичностью ≥ 0.90 (порог ниже, чем для A/B, так как
    запрос взят из A/B Langdon), отдельный столбец;
  * все локусы CBF сборки (16_04_cbf_loci.call_loci по окнам) — s16/copies/<g>_cbf_loci.tsv; локусы без копии-запроса Langdon
    перечисляются в столбце extra_cbf_loci (паралог по лучшему гену CS).
Выход: tables/16_gene_copies_by_genome.csv; genes16/<ген>_genomic.fa, _cds.fa, _protein.fa; s16/copies/all_copies.tsv."""
import os, re, sys, importlib.util
import numpy as np, pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("cl", os.path.join(HERE, "16_04_cbf_loci.py")); cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
spec = importlib.util.spec_from_file_location("cp", os.path.join(HERE, "16_04f_copies.py")); cp = importlib.util.module_from_spec(spec); spec.loader.exec_module(cp)
S16 = cl.WSL
GENOMES = ["Langdon", "Svevo", "Kronos", "ASM5281831v1", "Chili", "Mahmoudi", "Zavitan_WEW2.1", "Dicoccum_lpTriDico3",
           "CS_v1.0", "Jagger", "Julius", "SY_Mattis", "Norin61", "ArinaLrFor", "Stanley", "Landmark"]
LABEL = {"ASM5281831v1": "PI 192051", "Zavitan_WEW2.1": "Zavitan", "Dicoccum_lpTriDico3": "dicoccum lpTriDico3", "CS_v1.0": "Chinese Spring",
         "Stanley": "CDC Stanley", "Landmark": "CDC Landmark", "SY_Mattis": "SY Mattis", "Norin61": "Norin 61"}
SPECIES = {**{g: "твёрдая" for g in GENOMES[:6]}, "Zavitan_WEW2.1": "дикий эммер", "Dicoccum_lpTriDico3": "культурный эммер",
           **{g: "мягкая" for g in GENOMES[8:]}}
SV_REGION = "участок ±20 кб вокруг гена Svevo (синтения не проверялась)"
UNPLACED = "неразмещённый контиг"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s): return s.translate(COMP)[::-1]


def read_fa(path):
    d, n = {}, None
    if not os.path.exists(path): return d
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"): n = line[1:].strip().split()[0]; d[n] = []
        elif n: d[n].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}


def cut_win(wsq, r):
    """ген + 2 кб + 1 кб копии в ориентации гена из окна второго прохода (если копия лежит в окне)"""
    for k, s in wsq.items():
        m = re.match(r"^(.+):(\d+)-(\d+)$", k)
        if m and m.group(1) == r.tchr and int(m.group(2)) <= r.tg_s and int(m.group(3)) >= r.tg_e:
            o = int(m.group(2))
            a, b = (int(r.tg_s) - 2000, int(r.tg_e) + 1000) if r.strand == "+" else (int(r.tg_s) - 1000, int(r.tg_e) + 2000)
            x = s[max(0, a - o):b - o + 1]
            return x if r.strand == "+" else rc(x)
    return ""


def parse_mp(path):
    """miniprot GFF + ##STA -> mRNA с CDS и белком; координаты — хромосомные."""
    if not os.path.exists(path): return pd.DataFrame()
    rows, sta, cds = [], None, {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("##STA"): sta = line.rstrip("\n").split("\t")[1]; continue
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        m = re.match(r"^(.+):(\d+)-(\d+)$", f[0]); win, off = f[0], int(m.group(2)) - 1
        if f[2] == "mRNA":
            t, ts, te = a["Target"].split()
            rows.append(dict(id=a["ID"], win=win, chr=m.group(1), start=int(f[3]) + off, end=int(f[4]) + off, strand=f[6], score=float(f[5]),
                             ident=float(a["Identity"]), target=t, tstart=int(ts), tend=int(te), fs=int(a.get("Frameshift", 0)),
                             stop=int(a.get("StopCodon", 0)), prot=sta))
            sta = None
        elif f[2] == "CDS":
            cds.setdefault(a["Parent"], []).append((int(f[3]), int(f[4])))
    d = pd.DataFrame(rows)
    if len(d): d["cds"] = d.id.map(cds)
    return d


def cds_seq(r, winseq):
    s = winseq.get(r.win)
    if s is None or not isinstance(r.cds, list): return ""
    parts = [s[a - 1:b] for a, b in sorted(r.cds)]
    x = "".join(parts)
    return rc(x) if r.strand == "-" else x


def main(genomes):
    names = pd.read_csv(S16 + "q/candidate_names.tsv", sep="\t")
    cand = pd.read_csv(S16 + "q/candidates.tsv", sep="\t")
    cand = cand[cand.is_gene].merge(names[["cand_id", "ген", "приоритет", "локус"]], on="cand_id", how="left")
    plen = cl.read_fa_len(S16 + "q/prot2.fa"); plen.update(cl.read_fa_len(S16 + "q/prot_long.fa"))
    allc, summary, allloci = [], [], []
    for g in genomes:
        hp = S16 + f"copies/{g}_hits.tsv"
        if not os.path.exists(hp): print("нет", g); continue
        H = pd.read_csv(hp, sep="\t")
        H["genomic_seq"] = ""
        xp = S16 + f"p2/{g}_extra.paf"          # дополнительные запросы по синтенным окнам (16_04l_local_pass2.sh)
        if os.path.exists(xp) and os.path.getsize(xp):
            _, X = cp.analyse(g, xp, pd.read_csv(S16 + "q/candidates.tsv", sep="\t"))
            if len(X):
                Wn = pd.read_csv(S16 + f"copies/{g}_windows.tsv", sep="\t")
                gr = cand.set_index("cand_id").group
                X["synteny"] = [("в синтенном окне" if ((Wn.group == gr.get(r.cand_id)) & (Wn.tchr == r.tchr) & (Wn.wstart <= r.tg_s) &
                                                       (Wn.wend >= r.tg_e)).any() else "вне синтении") for _, r in X.iterrows()]
                wsq = read_fa(S16 + f"p2/{g}_win.fa")
                X["genomic_seq"] = [cut_win(wsq, r) for _, r in X.iterrows()]
                H = pd.concat([H, X], ignore_index=True)
        up = S16 + f"scan/{g}_un.paf"            # неразмещённые контиги (16_04n_server_unplaced.sh): синтения не оценивается
        if os.path.exists(up) and os.path.getsize(up):
            _, U = cp.analyse(g, up, pd.read_csv(S16 + "q/candidates.tsv", sep="\t"))
            if len(U):
                U["synteny"] = UNPLACED; U["genomic_seq"] = ""
                H = pd.concat([H, U], ignore_index=True)
        MP = pd.concat([parse_mp(S16 + f"p2/{g}_win.gff").assign(set="cbf"), parse_mp(S16 + f"p2/{g}_win_long.gff").assign(set="long")])
        if len(MP): MP["tlen"] = MP.target.map(plen); MP["cov"] = (MP.tend - MP.tstart + 1) / MP.tlen
        wins = read_fa(S16 + f"p2/{g}_win.fa"); cps = read_fa(S16 + f"p2/{g}_copies.fa")
        # все локусы CBF сборки
        if len(MP):
            d = cl.parse_gff(S16 + f"p2/{g}_win.gff", plen, chr_map=lambda s: s)
            loci = cl.call_loci(d, g) if len(d) else pd.DataFrame()
        else:
            loci = pd.DataFrame()
        loci.to_csv(S16 + f"copies/{g}_cbf_loci.tsv", sep="\t", index=False)
        for _, r in H.iterrows():
            c = cand[cand.cand_id == r.cand_id]
            if c.empty: continue
            c = c.iloc[0]
            is_cbf = r.cand_id.startswith("CBF_")
            m = MP[(MP.chr == r.tchr)] if len(MP) else MP
            if len(m):
                ov = (np.minimum(m.end, r.tg_e) - np.maximum(m.start, r.tg_s)).clip(lower=0) / (r.tg_e - r.tg_s + 1)
                m = m[ov >= 0.5]
                if is_cbf: m = m[m.set == "cbf"]
                else:
                    pref = m[m.target == f"SV|{r.cand_id}"]
                    m = pref if len(pref) else m[m.set == "long"]
            b = m.loc[m.score.idxmax()] if len(m) else None
            prot = (b.prot or "").rstrip("*") if b is not None else ""
            cds = cds_seq(b, wins) if b is not None else ""
            orf_ok = b is not None and b["cov"] >= 0.9 and b.fs == 0 and b.stop == 0
            why = [] if b is not None else ["нет выравнивания белка"]
            if b is not None:
                if b["cov"] < 0.9: why.append(f"покрытие белка {b['cov']:.2f}")
                if b.fs: why.append(f"сдвиг рамки ×{b.fs}")
                if b.stop: why.append(f"стоп ×{b.stop}")
            key = f"{r.cand_id}|{g}|{r.tchr}:{int(r.tg_s)}|{r.strand}"
            mchr = re.search(r"([1-7])([ABD])", str(r.tchr))
            if g in ("Chili", "Mahmoudi") or r.synteny == UNPLACED:   # скаффолды без надёжной хромосомы: гомеолог — по идентичности
                hom = c.chr[-1] if r.ident >= 0.95 else (mchr.group(2) if mchr else "?")
            else:
                hom = mchr.group(2) if mchr else "?"
            syn = r.synteny
            if g == "Svevo" and r.tchr in ("7B", "1B"):
                syn = SV_REGION
            allc.append(dict(genome=g, cand_id=r.cand_id, gene=c.ген, priority=c.приоритет, chr=r.tchr, start=int(r.tg_s), end=int(r.tg_e),
                             strand=r.strand, dna_ident=r.ident, gene_cov=r.gcov, synteny=syn, dna_complete=r.completeness,
                             homeolog=hom, indels=r.indels if isinstance(r.indels, str) else "",
                             mp_target=b.target if b is not None else "", mp_ident=b.ident if b is not None else None,
                             mp_cov=round(b["cov"], 3) if b is not None else None, orf="полная" if orf_ok else "неполная",
                             orf_why="; ".join(why), prot=prot, cds=cds, genomic=cps.get(key, "") or (r.genomic_seq if isinstance(r.genomic_seq, str) and r.genomic_seq else cut_win(wins, r))))
        print(g, "копий", len(H), "локусов CBF", len(loci))
        # все локусы CBF сборки: совпадение с копией гена Langdon (перекрытие) или «нет в Langdon»
        if len(loci):
            hh = H.merge(cand[["cand_id", "ген"]], on="cand_id", how="left")
            def match(r):
                x = hh[(hh.tchr == r.chr) & (hh.tg_s <= r.end) & (hh.tg_e >= r.start) & (hh.ident >= 0.95)]
                return ", ".join(sorted(set(x.ген.dropna()))) if len(x) else ""
            loci["ген_Langdon"] = [match(r) for _, r in loci.iterrows()]
            loci["genome"] = g
            allloci.append(loci)
    A = pd.DataFrame(allc)
    A.to_csv(S16 + "copies/all_copies.tsv", sep="\t", index=False)
    if allloci:      # все локусы CBF по miniprot в окнах группы 5 (у Chili/Mahmoudi — скаффолды)
        LL = pd.concat(allloci)
        LL["сборка"] = LL.genome.map(lambda g: LABEL.get(g, g))
        LL = LL[["сборка", "chr", "start", "end", "strand", "ген_Langdon", "cs_paralog", "cs_ident", "tm_name", "tm_ident", "best_target", "best_ident",
                 "best_cov", "complete", "why"]].rename(columns={"chr": "хромосома", "start": "начало", "end": "конец", "strand": "цепь",
                 "ген_Langdon": "ген_Langdon_по_ДНК", "cs_paralog": "лучший_ген_CS_той_же_хромосомы", "tm_name": "лучший_CBF_T_monococcum",
                 "complete": "полнота_белка", "why": "причина"})
        LL.to_csv(os.path.join(OUT, "tables/16_cbf_loci_by_genome.csv"), index=False, encoding="utf-8-sig")
    return A, cand


def tables(A, cand, genomes):
    os.makedirs(os.path.join(OUT, "genes16"), exist_ok=True)
    rows, fb = [], []
    for _, c in cand.iterrows():
        for g in genomes:
            x = A[(A.cand_id == c.cand_id) & (A.genome == g)]
            own = x[(x.homeolog == c.chr[-1]) & (x.dna_ident >= 0.95)]
            syn = own[own.synteny.isin(["в синтенном окне", SV_REGION])]
            outside = own[~own.synteny.isin(["в синтенном окне", SV_REGION, UNPLACED])]
            unpl = own[own.synteny == UNPLACED]
            D = x[(x.homeolog == "D") & (x.dna_ident >= 0.90)]
            f = lambda y: "; ".join(f"{r.chr}:{r.start}-{r.end}({r.strand})" for _, r in y.iterrows())
            rows.append(dict(ген=c.ген, cand_id=c.cand_id, локус=c.локус, приоритет=c.приоритет, сборка=LABEL.get(g, g), вид=SPECIES.get(g, ""),
                             найден="да" if len(syn) else ("вне синтении" if len(outside) else ("только на неразмещённом контиге" if len(unpl) else "нет")),
                             копий_в_синтенном_окне=len(syn), из_них_полных_ДНК=int((syn.dna_complete == "полная").sum()),
                             из_них_ORF_полная=int((syn.orf == "полная").sum()), копий_вне_синтении=len(outside), копий_на_неразмещённых_контигах=len(unpl),
                             координаты=f(syn), координаты_вне_синтении=f(outside), координаты_неразмещённые=f(unpl),
                             идентичность=";".join(f"{v:.4f}" for v in syn.dna_ident), покрытие_тела_гена=";".join(f"{v:.2f}" for v in syn.gene_cov),
                             ORF=";".join(o + (f" ({w})" if isinstance(w, str) and w else "") for o, w in zip(syn.orf, syn.orf_why)),
                             инделы_ДНК_10пн=";".join(v if isinstance(v, str) else "" for v in syn.indels),
                             гомеолог_D_копий=len(D) if SPECIES.get(g) == "мягкая" else None, гомеолог_D_координаты=f(D) if SPECIES.get(g) == "мягкая" else ""))
        # fasta
        base = re.sub(r"[^A-Za-z0-9_.-]+", "_", c.ген.replace(" @", "_")).strip("_")
        if c.cand_id.startswith("CBF_5A_513"): base += "_dup513"
        fb.append((base, c.cand_id))
        for kind, col in [("genomic", "genomic"), ("cds", "cds"), ("protein", "prot")]:
            x = A[(A.cand_id == c.cand_id) & (A.homeolog == c.chr[-1]) & (A.dna_ident >= 0.95) & (A[col].fillna("").map(len) > 0)
                  & A.synteny.isin(["в синтенном окне", SV_REGION])]          # в fasta — только синтенные копии (ортологи)
            x = x.assign(o=x.genome.map({g: i for i, g in enumerate(genomes)})).sort_values(["o", "start"])
            with open(os.path.join(OUT, "genes16", f"{base}_{kind}.fa"), "w", encoding="utf-8", newline="\n") as fh:
                k = {}
                for _, r in x.iterrows():
                    k[r.genome] = k.get(r.genome, 0) + 1
                    tag = "" if r.synteny in ("в синтенном окне", SV_REGION) else "|вне_синтении"
                    fh.write(f">{LABEL.get(r.genome, r.genome).replace(' ', '_')}|{r.chr}:{r.start}|копия{k[r.genome]}{tag}\n{r[col]}\n")
    T = pd.DataFrame(rows)
    pd.DataFrame(fb, columns=["base", "cand_id"]).to_csv(S16 + "q/file_base.tsv", sep="\t", index=False, header=False, lineterminator="\n")
    T.to_csv(os.path.join(OUT, "tables/16_gene_copies_by_genome.csv"), index=False, encoding="utf-8-sig")
    return T


if __name__ == "__main__":
    gl = sys.argv[1:] or GENOMES
    A, cand = main(gl)
    T = tables(A, cand, gl)
    pd.set_option("display.width", 250)
    print(T.pivot_table(index="ген", columns="сборка", values="копий_в_синтенном_окне", aggfunc="first").to_string())
