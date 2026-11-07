#!/usr/bin/env bash
# Этап 9.1. GWAS наборов по одному (тяжёлые задачи — последовательно), с проверкой памяти перед каждым.
source "$(dirname "$0")/00_env.sh"
set +e   # ошибка одного набора не должна обрывать остальные
for s in ${SETS:-RU IT S1 S2}; do
  echo "== $s $(date +%T)"; free -g | head -2
  "$ENVBIN/Rscript" "$(dirname "$0")/09_gwas.R" $s > "$WORK/logs/09_gwas_$s.log" 2>&1
  echo "exit $? ; $(grep -c '^[A-Z0-9]* .*: n=' "$WORK/logs/09_gwas_$s.log") прогонов"
done
