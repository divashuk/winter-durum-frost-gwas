#!/usr/bin/env bash
# Раздел 17, шаг 4д (служебный). Перенос готовых PAF сканирования с сервера (scan/<id>.done) в WSL s17/scan, окна второго прохода
# (17_08_pass2_regions.py, запускается из Windows отдельно) — обратно на сервер (q/pass2_<id>.tsv). Режимы: pull (PAF) | push (окна -> сервер, запуск 17_08b) | p2 (результаты второго прохода -> WSL).
set -uo pipefail
S=/work/scratch/s17; R=/work/remote/s17; mkdir -p $S/scan
if [ "$1" = pull ]; then
  for d in $(ssh compute-host "ls $R/scan/*.done 2>/dev/null"); do id=$(basename $d .done); [ "$id" = ALL ] && continue
    [ -s $S/scan/$id.paf ] || scp -q compute-host:$R/scan/$id.paf compute-host:$R/scan/${id}_un.paf $S/scan/ 2>/dev/null || scp -q compute-host:$R/scan/$id.paf $S/scan/; done
  ls $S/scan | wc -l
elif [ "$1" = p2 ]; then
  mkdir -p $S/p2
  for d in $(ssh compute-host "ls $R/p2/*.done 2>/dev/null"); do id=$(basename $d .done); [ "$id" = ALL ] && continue
    [ -s $S/p2/${id}_win_long.gff ] || for x in win.fa win.gff win_long.gff copies.fa; do scp -q compute-host:$R/p2/${id}_$x $S/p2/; done; done
  ls $S/p2 | wc -l
else
  scp -q $S/pass2/pass2_*.tsv compute-host:$R/q/ && ssh compute-host "bash $R/17_08b_server_pass2.sh"
fi
