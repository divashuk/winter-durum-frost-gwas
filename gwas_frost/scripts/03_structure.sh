#!/usr/bin/env bash
# Этап 3 (часть plink): LD-прунинг + PCA (plink2), попарный r2 внутри хромосом для спада LD (plink1.9).
# Использование: bash 03_structure.sh stream1|stream2
source "$(dirname "$0")/00_env.sh"
S=$1
Q="$WORK/qc/$S"; D="$WORK/pca/$S"; L="$WORK/ld/$S"; mkdir -p "$D" "$L"
LOG="$WORK/logs/03_structure_${S}.log"; : > "$LOG"
cd "$D"
plink2 --bfile "$Q/final" --allow-extra-chr --chr-set 14 no-xy --indep-pairwise 50 5 0.5 --out prune >> "$LOG" 2>&1
plink2 --bfile "$Q/final" --allow-extra-chr --chr-set 14 no-xy --extract prune.prune.in --pca 20 --out pca >> "$LOG" 2>&1
echo "$S: SNP после LD-прунинга (r2<0.5, окно 50): $(wc -l < prune.prune.in)" | tee -a "$LOG"
cd "$L"
plink --bfile "$Q/final" --allow-extra-chr --chr-set 14 --r2 --ld-window 999999 --ld-window-kb 100000 --ld-window-r2 0 \
      --threads $THREADS --out r2 >> "$LOG" 2>&1
echo "$S: пар r2 внутри хромосом: $(($(wc -l < r2.ld)-1))" | tee -a "$LOG"
gzip -f r2.ld
