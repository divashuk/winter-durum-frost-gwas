# -*- coding: utf-8 -*-
"""Раздел 17, часть A, шаг 4. Выбор основного гена-кандидата на локус -> tables/17_candidate_genes_other_loci.csv.
Гены и баллы: интервалы 14.4 (7B I05, 1B I06, 2B 767 I01, 2B 776 I02, 4A I07; tables/14_genes_all_intervals.csv) и новые интервалы
17 (1A 565.9 I09, 1A 597 I10, 4B 28.7 I11, 5B 502.8 I12; scripts/17_03_annot_score.py) — баллы scripts/14_score_genes.py.
Правило (задание 17, часть A):
 1) ген по критерию (а) (холод / яровизация / фотопериод) с r² ближайшего SNP к ведущему ≥ 0.2 — основной; если таких несколько —
    с наибольшим баллом; при равном балле — с большим r², затем ближе к ведущему SNP;
 2) иначе — ген, внутри которого лежит ведущий SNP; если такого нет — ген с наибольшим баллом среди генов без (а)
    (при равном балле — ближайший к ведущему SNP);
 «второй кандидат» — если основной выбран по правилу 2, а в блоке есть ген (а) с r² < 0.2: лучший из них (балл, r², расстояние).
 1B: дополнительно HSF TRITD1Bv1G203840 (по заданию).
Для 7B ведущий SNP — 7B:348683836 (IT, интервал I05); узкое окно I08 входит в I05."""
import os, sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
G = "/work/scratch/genes17/"
A = pd.read_csv(G + "genes_all9.csv")
IV = pd.read_csv(os.path.join(OUT, "tables/17_gene_intervals.csv"))
A = A[A.interval_id != "I08"]
DONE16 = {"TRITD1Bv1G205280", "TRITD7Bv1G107490", "TRITD7Bv1G098240"}     # уже разобраны в разделе 16 (16.4–16.7)
rows = []


def rec(r, role, rule, why):
    return dict(локус=r.locus, интервал=r.interval_id, ведущий_SNP=r.lead, панель=r.lead_set, роль=role, ген=r.svevo_gene,
                ID_CS_RefSeq_v1_1=r.cs_refseq_v11 if isinstance(r.cs_refseq_v11, str) else "",
                позиция_Langdon=f"{r.chr}:{int(r.langdon_start)}-{int(r.langdon_end)} ({r.strand})", расстояние_до_ведущего_кб=r.dist_to_lead_kb,
                r2_ближайшего_SNP=r.r2_nearest_with_lead, балл=int(r.score), баллы_а_б_в=f"{int(r.pts_a)}/{int(r.pts_b)}/{int(r.pts_c)}",
                критерий_а=r.criterion_a if isinstance(r.criterion_a, str) else "", функция=r.function_description, правило=rule, обоснование=why,
                разобран_в_разделе_16="да" if r.svevo_gene in DONE16 else "нет")


for iid, d in A.groupby("interval_id", sort=False):
    a = d[d.pts_a > 0]
    near = a[a.r2_nearest_with_lead >= 0.2].sort_values(["score", "r2_nearest_with_lead", "dist_to_lead_kb"], ascending=[False, False, True])
    far = a[a.r2_nearest_with_lead < 0.2].sort_values(["score", "r2_nearest_with_lead", "dist_to_lead_kb"], ascending=[False, False, True])
    n_a, n_near = len(a), len(near)
    if n_near:
        r = near.iloc[0]
        alt = "; ".join(f"{x.svevo_gene} (балл {int(x.score)}, r² {x.r2_nearest_with_lead})" for _, x in near.iloc[1:].iterrows())
        rows.append(rec(r, "основной", "1: ген (а) в блоке r² ≥ 0.2",
                        f"генов (а) в интервале {n_a}, из них с r² ≥ 0.2 — {n_near}; выбран по баллу, затем r², затем расстоянию"
                        + (f"; другие гены (а) с r² ≥ 0.2: {alt}" if alt else "")))
    else:
        inside = d[d.dist_to_lead_kb == 0]
        if len(inside):
            r = inside.sort_values("score", ascending=False).iloc[0]
            why = "ведущий SNP внутри гена"
        else:
            r = d[d.pts_a == 0].sort_values(["score", "dist_to_lead_kb"], ascending=[False, True]).iloc[0]
            why = "ведущий SNP вне генов; наибольший балл среди генов без (а), при равенстве — ближайший"
        why += f"; генов (а) в интервале {n_a}, с r² ≥ 0.2 — 0"
        rows.append(rec(r, "основной", "2: нет гена (а) с r² ≥ 0.2", why))
        if len(far):
            f = far.iloc[0]
            rows.append(rec(f, "второй кандидат", "ген (а) с r² < 0.2",
                            f"лучший из {len(far)} генов (а) интервала с r² < 0.2 (балл, r², расстояние)"))
    if iid == "I06":
        h = d[d.svevo_gene == "TRITD1Bv1G203840"].iloc[0]
        rows.append(rec(h, "дополнительно (по заданию)", "HSF по заданию", "HSF (IPR000232), 3.4 Мб от ведущего SNP; в 14.5 — без баллов (а)"))
T = pd.DataFrame(rows)
order = ["7B", "1B", "2B 767", "2B 775", "4A", "1A 565", "1A 597", "4B", "5B"]
T["_o"] = [next(i for i, o in enumerate(order) if l.startswith(o)) for l in T.локус]
T = T.sort_values(["_o"], kind="stable").drop(columns="_o")
T.to_csv(os.path.join(OUT, "tables/17_candidate_genes_other_loci.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(T[["локус", "роль", "ген", "позиция_Langdon", "расстояние_до_ведущего_кб", "r2_ближайшего_SNP", "балл", "критерий_а", "функция"]].to_string())
