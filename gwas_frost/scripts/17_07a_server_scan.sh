#!/bin/bash
# Раздел 17, шаг 4в (сервер compute-host; повтор 16_04a + 16_04c + 16_04n для генов части A). Рабочая папка — /work/remote/s17 (новая);
# папка раздела 16 /work/remote/s16 только читается (хромосомы групп 1, 5, 7, неразмещённые последовательности, заголовки сборок).
#  1) хромосомы групп 2 и 4 (тетраплоиды 2A 2B 4A 4B; мягкая 2A 2B 2D 4A 4B 4D) -> chr/<id>_g24.fa.gz — тем же кодом, что 16_04a
#     (имя хромосомы из заголовка); Chili и Mahmoudi — по всей сборке (индекс раздела 16);
#  2) minimap2 -x asm20 -I 4G --split-prefix -c --secondary=yes -N 50 -p 0.1 запросов (q/queries.fa: гены + якоря) по хромосомам групп 1, 5, 7
#     (/work/remote/s16/chr/<id>.fa.gz) и 2, 4 (chr/<id>_g24.fa.gz) -> scan/<id>.paf (объединение);
#  3) неразмещённые последовательности (/work/remote/s16/chr/<id>_un.fa.gz) — только гены (q/queries_Q.fa) -> scan/<id>_un.paf.
# Память: ожидание ≥ 20 ГБ доступной перед каждым minimap2, не больше 2 сборок одновременно (как в 16).
E=/opt/micromamba/envs/st_tools/bin
W=/work/remote/s17; W16=/work/remote/s16; mkdir -p $W/chr $W/scan $W/logs; cd $W
if [ "${1:-}" != run ]; then nohup setsid bash $W/17_07a_server_scan.sh run > logs/scan.log 2>&1 < /dev/null & echo started; exit 0; fi
waitmem() { while [ $(free -g | awk '/Mem/{print $7}') -lt 20 ]; do sleep 60; done; }
one() { id=$1; f=$2; [ -e scan/$id.done ] && return 0
  case $id in Langdon|ASM5281831v1|Kronos|Zavitan_WEW2.1|Dicoccum_lpTriDico3) T="2A 2B 4A 4B";; *) T="2A 2B 2D 4A 4B 4D";; esac
  if [ $id != Chili ] && [ $id != Mahmoudi ] && [ ! -e chr/$id.g24.done ]; then
    awk -F'\t' -v T="$T" 'BEGIN{n=split(T,a," "); for(i=1;i<=n;i++) w[a[i]]=1}
      { c=""; if (match($1, /chromosome:? ?[1-7][ABD]([, ]|$)/)) { s=substr($1,RSTART,RLENGTH); gsub(/chromosome:? ?|[, ]/,"",s); c=s }
        else { split($1,b," "); if (b[1] ~ /^[1-7][ABD]$/) c=b[1] }
        if (c in w) print $1 "\t" c }' $W16/headers/$id.tsv > chr/$id.names
    cut -f1 chr/$id.names > chr/$id.full
    case "$f" in
      *.tar.gz) tar -xzOf "$f" | $E/seqkit grep -n -f chr/$id.full ;;
      *.gz)     $E/seqkit grep -n -f chr/$id.full "$f" ;;
      *)        $E/samtools faidx "$f" $(awk '{print $1}' chr/$id.full) ;;
    esac | awk -F'\t' 'NR==FNR{m[$1]=$2; split($1,b," "); m[b[1]]=$2; next} /^>/{h=substr($0,2); print ">" m[h]; next} {print}' chr/$id.names - |
      $E/bgzip -@ 4 -c > chr/${id}_g24.fa.gz && $E/samtools faidx chr/${id}_g24.fa.gz && touch chr/$id.g24.done
  fi
  : > scan/$id.paf
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then TT="$f"; else TT="$W16/chr/$id.fa.gz chr/${id}_g24.fa.gz"; fi
  k=0; for t in $TT; do k=$((k+1)); waitmem
    $E/minimap2 -t 8 -x asm20 -I 4G --split-prefix scan/$id.tmp$k -c --secondary=yes -N 50 -p 0.1 "$t" q/queries.fa >> scan/$id.paf 2>> logs/mm_$id.log; done
  if [ -s $W16/chr/${id}_un.fa.gz ]; then waitmem
    $E/minimap2 -t 6 -x asm20 -c --secondary=yes -N 50 -p 0.1 $W16/chr/${id}_un.fa.gz q/queries_Q.fa > scan/${id}_un.paf 2> logs/mm_${id}_un.log; fi
  touch scan/$id.done; echo "$id $(date) $(wc -l < scan/$id.paf)"; }
export -f one waitmem; export E W W16
awk -F'\t' '$1!="CS_v2.1"' $W16/assemblies.tsv | xargs -P 2 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch scan/ALL.done
