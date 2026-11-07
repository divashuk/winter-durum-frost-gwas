# -*- coding: utf-8 -*-
"""Этап 10.3. Аннотация генов Svevo из широких интервалов (WSL genes10/genes_wide.csv) через Ensembl Plants BioMart
(датасет tturgidum_eg_gene): InterPro (ID + описание), GO (ID, название, домен), гомологи в T. aestivum (IWGSC, CS),
A. thaliana и O. sativa (Compara: тип гомологии, % идентичности, уверенность); описания генов Arabidopsis и риса
из athaliana_eg_gene / osativa_eg_gene. Сырые ответы — в WSL genes10/biomart_*.tsv. Собственных догадок нет."""
import io, os, sys, time, urllib.parse, urllib.request
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
W = r"/work/scratch/genes10"
URL = "https://plants.ensembl.org/biomart/martservice"
g = pd.read_csv(os.path.join(W, "genes_wide.csv"))
genes = sorted(g.svevo_gene.unique())
print("генов:", len(genes))


def mart(dataset, filt, ids, attrs):
    x = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="plants_mart" formatter="TSV" '
         'header="0" uniqueRows="1" datasetConfigVersion="0.6">'
         f'<Dataset name="{dataset}" interface="default"><Filter name="{filt}" value="{",".join(ids)}"/>'
         + "".join(f'<Attribute name="{a}"/>' for a in attrs) + "</Dataset></Query>")
    for k in range(5):
        try:
            with urllib.request.urlopen(URL + "?" + urllib.parse.urlencode({"query": x}), timeout=300) as r:   # POST сервер не принимает (405)
                txt = r.read().decode("utf-8")
            if "Query ERROR" in txt or "ERROR" in txt[:50]:
                raise RuntimeError(txt[:300])
            return pd.read_csv(io.StringIO(txt), sep="\t", header=None, names=attrs, dtype=str)
        except Exception as e:
            print("  повтор", k, str(e)[:200]); time.sleep(10 * (k + 1))
    raise RuntimeError(f"BioMart не ответил: {dataset} {attrs[:2]}")


def batched(dataset, filt, ids, attrs, name, n=60):
    f = os.path.join(W, f"biomart_{name}.tsv")
    old = pd.read_csv(f, sep="\t", dtype=str) if os.path.exists(f) else pd.DataFrame(columns=attrs)
    ids = [i for i in ids if i not in set(old.iloc[:, 0])]   # дозапрос только новых генов
    if not ids:
        return old
    parts = [old]
    for i in range(0, len(ids), n):
        parts.append(mart(dataset, filt, ids[i:i + n], attrs))
        print(f"  {name}: {min(i + n, len(ids))}/{len(ids)}")
    d = pd.concat(parts).drop_duplicates()
    d.to_csv(f, sep="\t", index=False)
    return d


ds, fl = "tturgidum_eg_gene", "ensembl_gene_id"
ipr = batched(ds, fl, genes, ["ensembl_gene_id", "interpro", "interpro_description"], "interpro")
go = batched(ds, fl, genes, ["ensembl_gene_id", "go_id", "name_1006", "namespace_1003"], "go")
hom = {}
for sp in ["taestivum", "athaliana", "osativa"]:
    a = ["ensembl_gene_id", f"{sp}_eg_homolog_ensembl_gene", f"{sp}_eg_homolog_associated_gene_name", f"{sp}_eg_homolog_orthology_type",
         f"{sp}_eg_homolog_perc_id", f"{sp}_eg_homolog_perc_id_r1", f"{sp}_eg_homolog_orthology_confidence", f"{sp}_eg_homolog_chromosome",
         f"{sp}_eg_homolog_chrom_start"]
    hom[sp] = batched(ds, fl, genes, a, f"hom_{sp}")
desc = {}
for sp in ["athaliana", "osativa", "taestivum"]:
    ids = sorted(hom[sp][f"{sp}_eg_homolog_ensembl_gene"].dropna().unique())
    desc[sp] = batched(f"{sp}_eg_gene", "ensembl_gene_id", ids, ["ensembl_gene_id", "external_gene_name", "description"], f"desc_{sp}", n=80)
for k, v in {"interpro": ipr, "go": go, **{f"hom_{s}": d for s, d in hom.items()}, **{f"desc_{s}": d for s, d in desc.items()}}.items():
    print(k, v.shape, "генов:", v.iloc[:, 0].nunique())
