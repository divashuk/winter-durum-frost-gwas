#!/usr/bin/env bash
# Раздел 17, часть C. Окружения для белков (WSL, GPU RTX 4070 12 ГБ; nvidia-smi: драйвер 610.53, CUDA 13.3). Окружение gwas_frost не меняется.
# Локальная сеть при установке обрывалась (pip: «Network is unreachable», conda: «< 30 байт/с»), поэтому пакеты скачивались на сервере compute-host
# и переносились scp: conda-пакеты — в кэш micromamba (~/micromamba/pkgs), wheels Python 3.12 — /work/remote/s17/wheels (srv_wheels.sh).
#  s17p  — Python 3.12, conda-forge/bioconda: TM-align, DSSP (mkdssp), PyMOL open-source, Biopython, pandas, matplotlib, scipy;
#          pip (офлайн, wheels): torch 2.7.1 (CUDA 12.6), transformers (ESMFold — EsmForProteinFolding; ESM-1v — EsmForMaskedLM), accelerate,
#          boltz[cuda] (Boltz-2);
#  cfold — Python 3.12 (venv на базе s17p): colabfold[alphafold] 1.6.1 + jax[cuda12] (версия по зависимостям colabfold); MSA — сервер MMseqs2 (api.colabfold.com).
# Веса моделей (сервер -> WSL, pull_weights.sh): ESM-1v 1–5 и ESMFold (Hugging Face, facebook/*) -> s17/models; AlphaFold2 (2021-07-14) ->
# ~/.cache/colabfold/params; Boltz-2 (boltz2_conf.ckpt, mols) -> ~/.boltz.   Лог: ~/oz_tv_sic_work/s17/logs/env_*.log
set -uo pipefail
E=/opt/micromamba/envs; L=/work/scratch/s17/logs; WH=/work/scratch/s17/wheels
mkdir -p $WH; scp -q -r compute-host:/work/remote/s17/wheels/torch compute-host:/work/remote/s17/wheels/cfold $WH/ 2>> $L/env_pip.log
$E/s17p/bin/pip install --no-index --find-links $WH/torch torch==2.7.1 transformers accelerate "boltz[cuda]" > $L/env_pip.log 2>&1
$E/s17p/bin/python -c "import torch, transformers, boltz; print('s17p', torch.__version__, torch.cuda.is_available(), transformers.__version__)" >> $L/env_pip.log 2>&1
[ -x $E/../cfold/bin/python ] || $E/s17p/bin/python -m venv /opt/micromamba/cfold
/opt/micromamba/cfold/bin/pip install --no-index --find-links $WH/cfold "colabfold[alphafold]==1.6.1" "jax[cuda12]" > $L/env_cfold.log 2>&1
/opt/micromamba/cfold/bin/python -c "import jax; print('cfold', jax.__version__, jax.devices())" >> $L/env_cfold.log 2>&1
tail -2 $L/env_pip.log $L/env_cfold.log
