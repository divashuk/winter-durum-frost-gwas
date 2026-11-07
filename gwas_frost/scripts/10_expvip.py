# -*- coding: utf-8 -*-
"""Этап 10.3. Холодовая экспрессия ортологов CS в expVIP (wheat-expression.com, набор генов RefSeq1.1).
Опыт SRP043554 (Li et al. 2015, Plant Cell, doi:10.1105/tpc.114.134338): побеги сорта Manitou, 3 листа,
холод 4 °C 2 недели против контроля, по 3 повтора. Для каждого гена: средний TPM холод / контроль,
log2((TPM_холод + 1) / (TPM_контроль + 1)), t-тест Уэлча. Ортолог CS — лучший по % идентичности Compara,
с приоритетом той же хромосомы. Один запрос возвращает всю триаду гомеологов. Кэш — WSL genes10/expvip_raw/."""
import json, os, sys, time, urllib.error, urllib.request
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
W = r"/work/scratch/genes10"
C = os.path.join(W, "expvip_raw"); os.makedirs(C, exist_ok=True)
STUDY = "SRP043554"
g = pd.read_csv(os.path.join(W, "genes_wide.csv"))
h = pd.read_csv(os.path.join(W, "biomart_hom_taestivum.tsv"), sep="\t", dtype=str).dropna(subset=["taestivum_eg_homolog_ensembl_gene"])
h["pid"] = h.taestivum_eg_homolog_perc_id.astype(float)
ch = g.drop_duplicates("svevo_gene").set_index("svevo_gene").svevo_chr
h["same_chr"] = [c == ch.get(s, "") for s, c in zip(h.ensembl_gene_id, h.taestivum_eg_homolog_chromosome)]
best = h.sort_values(["same_chr", "pid"], ascending=False).drop_duplicates("ensembl_gene_id")
best[["ensembl_gene_id", "taestivum_eg_homolog_ensembl_gene", "taestivum_eg_homolog_orthology_type", "pid", "same_chr",
      "taestivum_eg_homolog_chromosome", "taestivum_eg_homolog_chrom_start"]].to_csv(os.path.join(W, "cs_best_ortholog.tsv"), sep="\t", index=False)
# объём ограничен (сервер ~10 с на ген): сначала гены с функциональными критериями C1–C10, затем гены узких интервалов по удалённости от ведущего SNP
gw = pd.read_csv(r"/work/project/gwas_frost/tables/10_genes_wide.csv")
pri = pd.concat([gw[gw.func_score > 0].sort_values("dist_to_lead_kb"), gw[gw.in_narrow].sort_values("dist_to_lead_kb")]).svevo_gene.drop_duplicates()
m = best.set_index("ensembl_gene_id").taestivum_eg_homolog_ensembl_gene
todo = list(dict.fromkeys(m.reindex(pri).dropna()))
print("к запросу генов CS:", len(todo), flush=True)


def fetch(gid):
    f = os.path.join(C, gid + ".json")
    if (os.path.exists(f) and os.path.getsize(f) > 100) or os.path.exists(f + ".missing"):
        return gid, "cache"
    url = f"https://www.wheat-expression.com/expression_values/gene.json?compare=&gene_set=RefSeq1.1&name={gid}"
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                txt = r.read().decode("utf-8")
            d = json.loads(txt)
            groups = d["groups"]; exps = d["experiments"]
            keep = {e for e, x in exps.items() if groups[x["group"]]["factors"]["study"] == STUDY}
            out = {"gene": gid, "values": {}, "factor": {e: groups[exps[e]["group"]]["factors"]["Stress-disease"] for e in keep}}
            for tg, v in d["values"].items():
                out["values"][tg] = {e: x["value"] for e, x in v["tpm"].items() if e in keep}
            json.dump(out, open(f, "w"))
            return gid, "ok"
        except urllib.error.HTTPError as e:
            if e.code == 500:   # гена нет в наборе RefSeq1.1 expVIP — повтор не нужен
                open(f + ".missing", "w").close(); return gid, "нет в expVIP (HTTP 500)"
            err = str(e)[:100]; time.sleep(5 * (k + 1))
        except Exception as e:
            err = str(e)[:100]; time.sleep(5 * (k + 1))
    return gid, "ошибка: " + err


with ThreadPoolExecutor(3) as ex:
    for i, (gid, st) in enumerate(ex.map(fetch, todo)):
        if st.startswith("ошибка") or i % 200 == 0:
            print(i, gid, st, flush=True)
rows = []
for fn in [x for x in os.listdir(C) if x.endswith(".json")]:
    d = json.load(open(os.path.join(C, fn)))
    fac = d["factor"]
    for tg, v in d["values"].items():
        cold = [x for e, x in v.items() if fac[e] == "cold"]; ctl = [x for e, x in v.items() if fac[e] == "none"]
        if len(cold) and len(ctl):
            t = stats.ttest_ind(cold, ctl, equal_var=False) if len(cold) > 1 and len(ctl) > 1 else None
            rows.append({"cs_gene": tg, "tpm_cold": round(sum(cold) / len(cold), 2), "tpm_control": round(sum(ctl) / len(ctl), 2), "n_cold": len(cold), "n_control": len(ctl),
                         "p_welch": (t.pvalue if t is not None else None)})
E = pd.DataFrame(rows).drop_duplicates("cs_gene")
import numpy as np
E["log2fc_cold"] = np.round(np.log2((E.tpm_cold + 1) / (E.tpm_control + 1)), 2)
E.to_csv(os.path.join(W, "expvip_SRP043554.tsv"), sep="\t", index=False)
miss = [x[:-len(".json.missing")] for x in os.listdir(C) if x.endswith(".missing")]
pd.Series(miss, name="cs_gene").to_csv(os.path.join(W, "expvip_missing.tsv"), index=False)
print("нет в expVIP:", len(miss))
print("генов CS с данными:", E.shape[0], "; запрошено:", len(todo), "; найдено среди запрошенных:", E.cs_gene.isin(todo).sum())
