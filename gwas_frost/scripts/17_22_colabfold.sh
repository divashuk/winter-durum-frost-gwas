#!/usr/bin/env bash
# Раздел 17, часть C, шаг 2б (WSL, окружение cfold — venv /opt/micromamba/cfold, GPU). ColabFold (Mirdita et al. 2022; AlphaFold2, модель alphafold2_ptm, 1 модель (model_1; 5 моделей на 62 последовательности не укладывались по времени),
# 3 цикла рециркуляции — значения по умолчанию colabfold_batch) для аллелей выбранных пар: все гены CBF, PRR95 (TRITD5Bv1G171760)
# и гены, у пары которых есть индел > 5 а.о. (s17/prot/selected_pairs.tsv, столбец макс_индел_аа). MSA — сервер MMseqs2 (api.colabfold.com).
# Белки > 1000 а.о. не моделируются (как ESMFold). Выход: s17/prot/cf/out/<cand_id>_<аллель>_*rank_001*.pdb.
set -uo pipefail
E=/opt/micromamba/cfold/bin
P=/work/scratch/s17/prot; C=$P/cf; mkdir -p $C/out
$E/python - <<'PY'
import pandas as pd
P = "/work/scratch/s17/prot/"
SP = pd.read_csv(P + "selected_pairs.tsv", sep="\t")
seq, n = {}, None
for line in open(P + "alleles.fa"):
    if line.startswith(">"): n = line[1:].strip(); seq[n] = ""
    else: seq[n] += line.strip()
sel = SP[SP.cand_id.str.startswith("CBF_") | (SP.cand_id == "TRITD5Bv1G171760") | (SP.макс_индел_аа > 5)]
rows = [(f"{r.cand_id}_{al}", seq[f"{r.cand_id}|{al}"]) for _, r in sel.iterrows() for al in (r.a, r.b)]
rows = [(i, s) for i, s in rows if len(s) <= 1000]
with open(P + "cf/input.csv", "w") as f:
    f.write("id,sequence\n" + "".join(f"{i},{s}\n" for i, s in rows))
print("последовательностей", len(rows))
PY
export XLA_PYTHON_CLIENT_ALLOCATOR=platform TF_FORCE_GPU_ALLOW_GROWTH=true   # в WSL пинованная память хоста для JAX не выделялась (xla_gpu_host_bfc); с аллокатором platform — работает
$E/colabfold_batch --num-recycle 3 --num-models 1 --model-type alphafold2_ptm $C/a3m $C/out > $C/colabfold.log 2>&1
ls $C/out/*rank_001*.pdb | wc -l
