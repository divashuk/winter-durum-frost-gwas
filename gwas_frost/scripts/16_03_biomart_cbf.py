# -*- coding: utf-8 -*-
"""Раздел 16, шаг 3. Гены DREB1/CBF (InterPro IPR045277) в Svevo v1 (5A, 5B) и CS RefSeq v1.1 (5A, 5B, 5D) из Ensembl Plants BioMart;
гомологи Svevo -> CS (Compara) и InterPro всех генов-кандидатов. Только GET (POST сервер не принимает, см. 10_biomart.py).
Ответы сохраняются в WSL ~/oz_tv_sic_work/s16/biomart/*.tsv; собственных догадок нет."""
import io, os, sys, time, urllib.parse, urllib.request
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
W = "/work/scratch/s16/biomart"
os.makedirs(W, exist_ok=True)
URL = "https://plants.ensembl.org/biomart/martservice"


def mart(dataset, filters, attrs):
    f = "".join(f'<Filter name="{k}" value="{v}"/>' for k, v in filters.items())
    x = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="plants_mart" formatter="TSV" '
         'header="0" uniqueRows="1" datasetConfigVersion="0.6">'
         f'<Dataset name="{dataset}" interface="default">{f}' + "".join(f'<Attribute name="{a}"/>' for a in attrs) + "</Dataset></Query>")
    for k in range(5):
        try:
            with urllib.request.urlopen(URL + "?" + urllib.parse.urlencode({"query": x}), timeout=300) as r:
                txt = r.read().decode("utf-8")
            if "ERROR" in txt[:200]:
                raise RuntimeError(txt[:300])
            return pd.read_csv(io.StringIO(txt), sep="\t", header=None, names=attrs, dtype=str)
        except Exception as e:
            print("  повтор", k, str(e)[:200]); time.sleep(10 * (k + 1))
    raise RuntimeError("BioMart не ответил")


def save(d, name):
    d.to_csv(os.path.join(W, name + ".tsv"), sep="\t", index=False)
    print(name, len(d))
    return d


loc = ["ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"]
sv = save(mart("tturgidum_eg_gene", {"interpro": "IPR045277", "chromosome_name": "5A,5B"}, loc + ["description"]), "svevo_IPR045277_5A5B")
cs = save(mart("taestivum_eg_gene", {"interpro": "IPR045277", "chromosome_name": "5A,5B,5D"}, loc + ["external_gene_name", "description"]),
          "cs_IPR045277_5ABD")
# гомологи Svevo -> CS для генов DREB1 и прочих кандидатов
cand = sorted(set(sv.ensembl_gene_id) | {"TRITD5Av1G171440", "TRITD5Bv1G171760", "TRITD7Bv1G107490", "TRITD7Bv1G098240", "TRITD1Bv1G205280",
                                         "TRITD5Bv1G171430", "TRITD5Bv1G171580", "TRITD5Bv1G171610"})
hom = save(mart("tturgidum_eg_gene", {"ensembl_gene_id": ",".join(cand)},
                ["ensembl_gene_id", "taestivum_eg_homolog_ensembl_gene", "taestivum_eg_homolog_orthology_type", "taestivum_eg_homolog_perc_id",
                 "taestivum_eg_homolog_perc_id_r1", "taestivum_eg_homolog_chromosome", "taestivum_eg_homolog_chrom_start"]), "svevo_cs_homologs")
ipr = save(mart("tturgidum_eg_gene", {"ensembl_gene_id": ",".join(cand)}, ["ensembl_gene_id", "interpro", "interpro_description",
                                                                         "interpro_start", "interpro_end"]), "svevo_interpro")
pos = save(mart("tturgidum_eg_gene", {"ensembl_gene_id": ",".join(cand)}, loc + ["description"]), "svevo_positions")
csi = save(mart("taestivum_eg_gene", {"ensembl_gene_id": ",".join(sorted(cs.ensembl_gene_id))},
                ["ensembl_gene_id", "interpro", "interpro_description", "interpro_start", "interpro_end"]), "cs_interpro")
