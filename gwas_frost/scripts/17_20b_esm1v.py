# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 1б (WSL, окружение s17p, GPU). Оценка замен ESM-1v (Meier et al. 2021): масочная лог-вероятность
(masked marginal): позиция i последовательности аллеля-фона маскируется, score = log P(мутантный а.о.) − log P(исходный а.о.);
ансамбль 5 моделей esm1v_t33_650M_UR90S_1…5 (веса Hugging Face facebook/esm1v_t33_650M_UR90S_k — те же модели, что в fair-esm;
transformers EsmForMaskedLM), итог — среднее по моделям. Белки длиннее 1022 а.о. — окно 1022 а.о. с позицией в центре (у краёв — сдвинутое).
Задания — s17/prot/esm_jobs.tsv (17_20a): каждая замена пары в обе стороны (фон a и фон b). Выход: s17/prot/esm_scores.tsv."""
import time
import pandas as pd, torch
from transformers import AutoTokenizer, EsmForMaskedLM

P = "/work/scratch/s17/prot/"; MD = "/work/scratch/s17/models/"
J = pd.read_csv(P + "esm_jobs.tsv", sep="\t")
seq, n = {}, None
for line in open(P + "alleles.fa"):
    if line.startswith(">"): n = line[1:].strip(); seq[n] = ""
    else: seq[n] += line.strip()
dev = "cuda" if torch.cuda.is_available() else "cpu"
print("устройство", dev, "заданий", len(J), flush=True)
res = {}
for k in range(1, 6):
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(MD + f"esm1v_{k}")
    model = EsmForMaskedLM.from_pretrained(MD + f"esm1v_{k}").eval().to(dev)
    if dev == "cuda": model = model.half()
    for (cid, sid, pos), g in J.groupby(["cand_id", "seq_id", "pos"]):
        s = seq[f"{cid}|{sid}"]; i = int(pos) - 1
        if len(s) > 1022:
            lo = min(max(0, i - 511), len(s) - 1022); s = s[lo:lo + 1022]; i -= lo
        x = tok(s, return_tensors="pt")["input_ids"]
        x[0, i + 1] = tok.mask_token_id                       # +1: токен <cls> в начале
        with torch.no_grad():
            lp = torch.log_softmax(model(input_ids=x.to(dev)).logits[0, i + 1].float(), -1).cpu()
        for _, r in g.iterrows():
            assert s[i] == r.wt, (cid, sid, pos, s[i], r.wt)
            res.setdefault((cid, sid, int(pos), r.mt), []).append(float(lp[tok.convert_tokens_to_ids(r.mt)] - lp[tok.convert_tokens_to_ids(r.wt)]))
    print("модель", k, round(time.time() - t0), "с", flush=True)
    del model; torch.cuda.empty_cache()
out = [dict(cand_id=c, seq_id=s, pos=p, mt=m, **{f"esm1v_{i + 1}": v[i] for i in range(5)}, esm1v_mean=sum(v) / len(v)) for (c, s, p, m), v in res.items()]
J.merge(pd.DataFrame(out), on=["cand_id", "seq_id", "pos", "mt"]).to_csv(P + "esm_scores.tsv", sep="\t", index=False)
print("готово", len(out))
