# -*- coding: utf-8 -*-
"""Этап 5.4. Функциональная аннотация генов Svevo в интервалах QTL — только из Ensembl Plants BioMart
(датасет tturgidum_eg_gene): InterPro-домены и GO-термины. В GFF3/кДНК Ensembl описаний генов Svevo нет.
Собственные догадки о функции не добавляются. Дополняет tables/05_genes_in_qtl.csv столбцами interpro, go_terms."""
import os, sys, time, urllib.parse, urllib.request, io
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
T = r"/work/project/gwas_frost/tables"
g = pd.read_csv(os.path.join(T, "05_genes_in_qtl.csv"))
genes = sorted(g.svevo_gene.unique())
URL = "https://plants.ensembl.org/biomart/martservice"

def mart(ids, attr):
    x = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="plants_mart" formatter="TSV" '
         'header="0" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="tturgidum_eg_gene" interface="default">'
         f'<Filter name="ensembl_gene_id" value="{",".join(ids)}"/><Attribute name="ensembl_gene_id"/>'
         f'<Attribute name="{attr}"/></Dataset></Query>')
    for k in range(4):
        try:
            with urllib.request.urlopen(URL + "?" + urllib.parse.urlencode({"query": x}), timeout=180) as r:
                txt = r.read().decode("utf-8")
            if "Query ERROR" in txt:
                raise RuntimeError(txt[:200])
            return pd.read_csv(io.StringIO(txt), sep="\t", header=None, names=["gene", "term"], dtype=str)
        except Exception as e:
            print("повтор", k, e); time.sleep(5)
    return pd.DataFrame(columns=["gene", "term"])

out = {}
for attr, col in [("interpro_description", "interpro"), ("name_1006", "go_terms")]:
    parts = [mart(genes[i:i + 60], attr) for i in range(0, len(genes), 60)]
    d = pd.concat(parts).dropna()
    out[col] = d.groupby("gene").term.apply(lambda s: "; ".join(sorted(set(s))))
    print(attr, "генов с термами:", out[col].shape[0], "из", len(genes))
g = g.drop(columns=[c for c in ["interpro", "go_terms"] if c in g.columns])
g["interpro"] = g.svevo_gene.map(out["interpro"]).fillna("")
g["go_terms"] = g.svevo_gene.map(out["go_terms"]).fillna("")
g["annotation"] = g.apply(lambda r: r.interpro if r.interpro else ("нет InterPro/GO в Ensembl Plants" if not r.go_terms else "только GO"), axis=1)
g.to_csv(os.path.join(T, "05_genes_in_qtl.csv"), index=False, encoding="utf-8-sig")
print(g.groupby("QTL").agg(n=("svevo_gene", "size"), with_interpro=("interpro", lambda s: (s != "").sum())).to_string())
