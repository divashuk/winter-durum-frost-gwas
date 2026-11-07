#!/bin/bash
# Раздел 16, шаг 4о (сервер compute-host). Неразмещённые последовательности сборок хромосомного уровня (всё, кроме псевдомолекул групп
# 1–7; Chili и Mahmoudi уже просмотрены целиком): извлечение -> chr/<id>_un.fa.gz; minimap2 -x asm20 только генов-кандидатов
# (q/queries_Q.fa, без якорей) -> scan/<id>_un.paf. Копии на них — «неразмещённые контиги» (синтения не оценивается).
E=/opt/micromamba/envs/st_tools/bin
W=/work/remote/s16; cd $W
if [ "${1:-}" != run ]; then nohup setsid bash $W/16_04n_server_unplaced.sh run > logs/unplaced.log 2>&1 < /dev/null & echo started; exit 0; fi
while [ ! -e scan/ALL.done ]; do sleep 60; done
waitmem() { while [ $(free -g | awk '/Mem/{print $7}') -lt 20 ]; do sleep 60; done; }
one() { id=$1; f=$2; [ -e scan/${id}_un.done ] && return 0
  case $id in Chili|Mahmoudi|CS_v2.1) return 0;; esac
  awk -F'\t' '{ if ($1 ~ /chromosome:? ?[1-7][ABD]([, ]|$)/ || $1 ~ /^[1-7][ABD]( |$)/) next; print $1 }' headers/$id.tsv > chr/$id.unnames
  case "$f" in
    *.tar.gz) tar -xzOf "$f" | $E/seqkit grep -n -f chr/$id.unnames ;;
    *.gz)     $E/seqkit grep -n -f chr/$id.unnames "$f" ;;
    *)        $E/samtools faidx "$f" $(awk '{print $1}' chr/$id.unnames) ;;
  esac | $E/seqkit replace -p '^(\S+).*' -r '$1' | $E/bgzip -@ 2 -c > chr/${id}_un.fa.gz
  waitmem
  $E/minimap2 -t 6 -x asm20 -c --secondary=yes -N 50 -p 0.1 chr/${id}_un.fa.gz q/queries_Q.fa > scan/${id}_un.paf 2> logs/mm_${id}_un.log && touch scan/${id}_un.done; }
export -f one waitmem; export E
awk -F'\t' '{print $1 "\t" $2}' assemblies.tsv | xargs -P 3 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch scan/UN.done
