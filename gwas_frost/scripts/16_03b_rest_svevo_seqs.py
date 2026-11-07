# -*- coding: utf-8 -*-
"""Раздел 16, шаг 3б. Белок, CDS и кДНК канонического транскрипта генов-кандидатов Svevo v1 (Ensembl REST, rest.ensembl.org).
Кандидаты: все IPR045277 Svevo на 5A/5B (biomart/svevo_IPR045277_5A5B.tsv) + PRR95, TRITD5Av1G171440, 7B, 1B и гены у SNP 5B.
Выход: WSL ~/oz_tv_sic_work/s16/seq/svevo_{pep,cds,cdna}.fa и svevo_canonical.tsv."""
import json, os, sys, time, urllib.request
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
B = "/work/scratch/s16/"
os.makedirs(B + "seq", exist_ok=True)
R = "https://rest.ensembl.org"


def get(path, ctype="application/json"):
    for k in range(6):
        try:
            req = urllib.request.Request(R + path, headers={"Content-Type": ctype})
            with urllib.request.urlopen(req, timeout=120) as r:
                t = r.read().decode()
            time.sleep(0.07)
            return json.loads(t) if ctype == "application/json" else t
        except Exception as e:
            print("  повтор", k, path, str(e)[:120]); time.sleep(3 * (k + 1))
    raise RuntimeError(path)


sv = pd.read_csv(B + "biomart/svevo_IPR045277_5A5B.tsv", sep="\t")
genes = list(sv.ensembl_gene_id) + ["TRITD5Av1G171440", "TRITD5Bv1G171760", "TRITD7Bv1G107490", "TRITD7Bv1G098240",
                                    "TRITD1Bv1G205280", "TRITD5Bv1G171580", "TRITD5Bv1G171610", "TRITD5Bv1G171430"]
from concurrent.futures import ThreadPoolExecutor


def one(g):   # запросы одного гена; гены — параллельно (сервер отвечает ~10–15 с на запрос)
    j = get(f"/lookup/id/{g}?expand=1")
    tr = [t for t in j["Transcript"] if t.get("is_canonical")] or j["Transcript"]
    t = tr[0]["id"]
    row = dict(gene=g, transcript=t, chr=j["seq_region_name"], start=j["start"], end=j["end"], strand=j["strand"],
               n_exons=len(tr[0].get("Exon", [])))
    seqs = {typ: get(f"/sequence/id/{t}?type={typ}", "text/plain") for typ in ["protein", "cds", "cdna"]}
    print(g, t, flush=True)
    return row, seqs


rows, fa = [], {"pep": [], "cds": [], "cdna": []}
with ThreadPoolExecutor(6) as ex:
    for row, seqs in ex.map(one, genes):
        rows.append(row)
        for typ, s in seqs.items():
            fa["pep" if typ == "protein" else typ].append(f">{row['gene']} {row['transcript']}\n{s}\n")
for k, v in fa.items():
    open(B + f"seq/svevo_{k}.fa", "w", newline="\n").write("".join(v))
pd.DataFrame(rows).to_csv(B + "seq/svevo_canonical.tsv", sep="\t", index=False, lineterminator="\n")
