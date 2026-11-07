# -*- coding: utf-8 -*-
"""Раздел 16, шаг 6а. Ортологи генов-кандидатов в ячмене (Hordeum vulgare), рисе (Oryza sativa) и Brachypodium distachyon
по Ensembl Compara (rest.ensembl.org /homology/id, compara=plants, type=orthologues, sequence=protein).
Для гена-кандидата берётся ID Svevo (если ген есть в Svevo: лучший белок Svevo с идентичностью ≥ 0.97), иначе — лучший ген CS (IWGSC v1.1).
Белки ортологов сохраняются (без выравнивания: выравнивание по позициям делается MAFFT в шаге 5/6).
Выход: WSL s16/compara/orthologs.tsv, orthologs_pep.fa."""
import json, os, sys, time, urllib.request
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

sys.stdout.reconfigure(encoding="utf-8")
S16 = "/work/scratch/s16/"
os.makedirs(S16 + "compara", exist_ok=True)
R = "https://rest.ensembl.org"
TARGET = ["hordeum_vulgare", "oryza_sativa", "brachypodium_distachyon"]


def get(path):
    for k in range(6):
        try:
            req = urllib.request.Request(R + path, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print("  повтор", k, path[:80], str(e)[:100], flush=True); time.sleep(4 * (k + 1))
    return None


L = pd.read_csv(S16 + "langdon_cbf/langdon_cbf_loci.tsv", sep="\t")
C = pd.read_csv(S16 + "q/candidates.tsv", sep="\t")
ids = {}
for _, r in C[C.is_gene].iterrows():
    if r.cand_id.startswith("TRITD"):
        ids[r.cand_id] = ("triticum_turgidum", r.cand_id)
    else:
        l = L[(L.chr == r.chr) & (L.start == r.start)].iloc[0]
        if isinstance(l.svevo_best, str) and l.svevo_ident >= 0.97:
            ids[r.cand_id] = ("triticum_turgidum", l.svevo_best)
        else:
            ids[r.cand_id] = ("triticum_aestivum", l.cs_paralog)


def one(item):
    cand, (sp, gid) = item
    q = "&".join(f"target_species={t}" for t in TARGET)
    j = get(f"/homology/id/{sp}/{gid}?compara=plants;type=orthologues;sequence=protein;aligned=0;{q.replace('&', ';')}")
    rows, fa = [], []
    for h in (j or {}).get("data", [{}])[0].get("homologies", []) if j else []:
        t = h["target"]
        rows.append(dict(cand_id=cand, query_species=sp, query_gene=gid, target_species=t["species"], target_gene=t["id"],
                         target_protein=t["protein_id"], type=h["type"], perc_id=t["perc_id"], perc_pos=t["perc_pos"]))
        fa.append(f">{cand}|{t['species']}|{t['id']}\n{t['seq']}\n")
    print(cand, gid, len(rows), flush=True)
    return rows, fa


rows, fa = [], []
with ThreadPoolExecutor(6) as ex:
    for r_, f_ in ex.map(one, ids.items()):
        rows += r_; fa += f_
pd.DataFrame(rows).to_csv(S16 + "compara/orthologs.tsv", sep="\t", index=False, lineterminator="\n")
open(S16 + "compara/orthologs_pep.fa", "w", newline="\n").write("".join(fa))
pd.DataFrame([dict(cand_id=k, species=v[0], gene=v[1]) for k, v in ids.items()]).to_csv(S16 + "compara/query_ids.tsv", sep="\t", index=False,
                                                                                       lineterminator="\n")
