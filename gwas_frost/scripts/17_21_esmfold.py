# -*- coding: utf-8 -*-
"""Раздел 17, часть C, шаг 2а (WSL, окружение s17p, GPU RTX 4070 12 ГБ). 3D-модели обоих аллелей выбранной пары (s17/prot/selected_pairs.tsv)
ESMFold (Lin et al. 2023; transformers EsmForProteinFolding, веса facebook/esmfold_v1 — s17/models/esmfold_v1) для белков ≤ 1000 а.о.
Ствол ESM-2 — в half precision, свёртка — float32 (белки ≤ 400 а.о.) или bfloat16 (> 400 а.о., запуск с аргументом long:
иначе модель не помещалась в свободную память GPU); chunk_size по длине белка: ≤ 400 а.о. — 64, ≤ 700 — 32, длиннее — 16 (при нехватке памяти GPU — меньше; WSL при нехватке
памяти GPU уходит в общую память и замедляется в десятки раз). pLDDT пишется в B-фактор (0–100).
Выход: s17/prot/esmfold/<cand_id>_<аллель>.pdb, esmfold_log.tsv."""
import os, sys, time
import pandas as pd, torch
from transformers import AutoTokenizer, EsmForProteinFolding

P = "/work/scratch/s17/prot/"; O = P + "esmfold/"; os.makedirs(O, exist_ok=True)
SP = pd.read_csv(P + "selected_pairs.tsv", sep="\t")
seq, n = {}, None
for line in open(P + "alleles.fa"):
    if line.startswith(">"): n = line[1:].strip(); seq[n] = ""
    else: seq[n] += line.strip()
tok = AutoTokenizer.from_pretrained("/work/scratch/s17/models/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("/work/scratch/s17/models/esmfold_v1", low_cpu_mem_usage=True)
LONG = len(sys.argv) > 1 and sys.argv[1] == "long"      # второй проход: белки > 400 а.о., блок свёртки в bfloat16 (иначе не помещается в память GPU)
if LONG: model = model.to(torch.bfloat16)
model.esm = model.esm.half(); model = model.cuda().eval()
torch.backends.cuda.matmul.allow_tf32 = True
log = []
jobs = [(r.cand_id, x) for _, r in SP.iterrows() for x in (r.a, r.b)]
for cid, al in jobs:
    s = seq[f"{cid}|{al}"]; out = O + f"{cid}_{al}.pdb"
    if len(s) > 1000: log.append(dict(cand_id=cid, allele=al, length=len(s), status="> 1000 а.о. — ESMFold не запускался")); continue
    if os.path.exists(out): log.append(dict(cand_id=cid, allele=al, length=len(s), status="готово ранее")); continue
    if (len(s) > 400) != LONG: continue
    t0 = time.time(); st = "ошибка"
    for ch in [c for c in (64, 32, 16, 8) if c <= (64 if len(s) <= 400 else 32 if len(s) <= 700 else 16)]:
        try:
            model.trunk.set_chunk_size(ch)
            x = tok([s], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
            with torch.no_grad(): o = model(x)
            o = {k: (v.float() if torch.is_tensor(v) and v.is_floating_point() else v) for k, v in o.items()}
            pdb = model.output_to_pdb(o)[0]
            lines = []
            for l in pdb.splitlines():          # B-фактор: pLDDT 0–1 -> 0–100
                if l.startswith("ATOM") and float(l[60:66]) <= 1.0: l = l[:60] + f"{float(l[60:66]) * 100:6.2f}" + l[66:]
                lines.append(l)
            open(out, "w").write("\n".join(lines) + "\n")
            st = f"ok ({'bf16, ' if LONG else ''}chunk {ch}, {time.time() - t0:.0f} с, pLDDT {float(o['plddt'][0, :, 1].float().mean()) * 100:.1f})"; break
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache(); st = f"нехватка памяти GPU при chunk {ch}"
    log.append(dict(cand_id=cid, allele=al, length=len(s), status=st)); print(cid, al, len(s), st, flush=True)
    torch.cuda.empty_cache()
pd.DataFrame(log).to_csv(P + "esmfold_log.tsv", sep="\t", index=False)
