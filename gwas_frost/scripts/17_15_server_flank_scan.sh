#!/bin/bash
# Раздел 17, шаг 7ж (сервер compute-host). Аллели SNP типа (б) (разделы 16 и 17) в 16 сборках: шаблоны Langdon ± 150 п.н. (q/flanks_b.fa)
# -> minimap2 -x asm20 -c --cs=long --secondary=yes -N 20 по хромосомам групп 1, 2, 4, 5, 7 (Chili, Mahmoudi — вся сборка; без --split-prefix: он несовместим с --cs) и неразмещённым
# последовательностям -> scan_b/<id>.paf. Аллель сборки = основание мишени в позиции 151 шаблона (разбор — 17_15b_flank_alleles.py).
E=/opt/micromamba/envs/st_tools/bin
W=/work/remote/s17; W16=/work/remote/s16; mkdir -p $W/scan_b $W/logs; cd $W
if [ "${1:-}" != run ]; then nohup setsid bash $W/17_15_server_flank_scan.sh run > logs/scan_b.log 2>&1 < /dev/null & echo started; exit 0; fi
waitmem() { while [ $(free -g | awk '/Mem/{print $7}') -lt 20 ]; do sleep 60; done; }
one() { id=$1; f=$2; [ -e scan_b/$id.done ] && return 0
  while [ ! -e scan/$id.done ]; do sleep 60; done          # хромосомы групп 2, 4 готовит 17_07a
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then TT="$f"; else TT="$W16/chr/$id.fa.gz chr/${id}_g24.fa.gz $W16/chr/${id}_un.fa.gz"; fi
  : > scan_b/$id.paf; k=0
  for t in $TT; do [ -s "$t" ] || continue; k=$((k+1)); waitmem
    $E/minimap2 -t 6 -x asm20 -I 4G -c --cs=long --secondary=yes -N 20 -p 0.5 "$t" q/flanks_b.fa >> scan_b/$id.paf 2>> logs/mmb_$id.log; done
  touch scan_b/$id.done; echo "$id $(date)"; }
export -f one waitmem; export E W W16
awk -F'\t' '$1!="CS_v2.1"' $W16/assemblies.tsv | xargs -P 3 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch scan_b/ALL.done
