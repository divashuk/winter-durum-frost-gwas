#!/usr/bin/env bash
# Этап 4. Запуск GWAS обоих потоков внутри WSL (отвязано от клиента wsl.exe); ждёт завершения обоих.
source "$(dirname "$0")/00_env.sh"
SC="$(dirname "$0")/04_gwas.R"
for s in stream1 stream2; do
  setsid nohup "$ENVBIN/Rscript" "$SC" $s > "$WORK/logs/04_gwas_$s.log" 2>&1 &
done
wait
grep -h "^stream\|rror" "$WORK"/logs/04_gwas_stream*.log | tail -60
