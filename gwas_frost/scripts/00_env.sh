#!/usr/bin/env bash
# Общие пути и окружение для всех этапов (source этот файл в начале каждого скрипта).
# Тяжёлые файлы — только в WORK (вне OneDrive).
set -euo pipefail
export ENVBIN=/opt/micromamba/envs/gwas_frost/bin
export PATH="$ENVBIN:$PATH"
export WORK=/work/scratch
export PROJ="/work/project"
export OUT="$PROJ/gwas_frost"
export SRC="/work/phenotypes"
export THREADS=12
mkdir -p "$WORK"/{raw,qc,imp,pca,kin,ld,gwas,logs,tmp}
