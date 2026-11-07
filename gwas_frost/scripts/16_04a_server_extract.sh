#!/bin/bash
# Раздел 16, шаг 4а (сервер compute-host). Целевые хромосомы сборок -> /work/remote/s16/chr/<id>.fa.gz (bgzip + .fai/.gzi).
# Хранилище /storage1 только читается. Тетраплоиды: 1A 1B 5A 5B 7A 7B; мягкая: 1A 1B 1D 5A 5B 5D 7A 7B 7D.
# Chili и Mahmoudi (RagTag, хромосомы — тысячи скаффолдов, приписка к хромосомам ненадёжна) не копируются:
# строится только индекс BGZF (.fai/.gzi) в /work/remote/s16/chr/, поиск — по всей сборке.
# Имя хромосомы — из заголовка: «chromosome: 5A», «chromosome 5A,», или первое слово «5A».
# Запуск: ssh compute-host 'bash -s' < 16_04a_server_extract.sh   (уходит в фон, флаги chr/<id>.done)
E=/opt/micromamba/envs/st_tools/bin
W=/work/remote/s16; mkdir -p $W/chr $W/logs; cd $W
if [ "${1:-}" != run ]; then nohup setsid bash $W/16_04a_server_extract.sh run > logs/extract.log 2>&1 < /dev/null & echo started; exit 0; fi
while [ ! -e headers/ALL.done ]; do sleep 30; done
one() { id=$1; f=$2; [ -e chr/$id.done ] && return 0
  case $id in Langdon|ASM5281831v1|Kronos|Zavitan_WEW2.1|Dicoccum_lpTriDico3) T="1A 1B 5A 5B 7A 7B";; *) T="1A 1B 1D 5A 5B 5D 7A 7B 7D";; esac
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then
    $E/samtools faidx --fai-idx chr/$id.fai --gzi-idx chr/$id.gzi "$f" && touch chr/$id.done; return; fi
  awk -F'\t' -v T="$T" 'BEGIN{n=split(T,a," "); for(i=1;i<=n;i++) w[a[i]]=1}
    { c=""; if (match($1, /chromosome:? ?[1-7][ABD]([, ]|$)/)) { s=substr($1,RSTART,RLENGTH); gsub(/chromosome:? ?|[, ]/,"",s); c=s }
      else { split($1,b," "); if (b[1] ~ /^[1-7][ABD]$/) c=b[1] }
      if (c in w) print $1 "\t" c }' headers/$id.tsv > chr/$id.names
  cut -f1 chr/$id.names > chr/$id.full
  case "$f" in
    *.tar.gz) tar -xzOf "$f" | $E/seqkit grep -n -f chr/$id.full ;;
    *.gz)     $E/seqkit grep -n -f chr/$id.full "$f" ;;
    *)        $E/samtools faidx "$f" $(awk '{print $1}' chr/$id.full) ;;
  esac | awk -F'\t' 'NR==FNR{m[$1]=$2; split($1,b," "); m[b[1]]=$2; next} /^>/{h=substr($0,2); print ">" m[h]; next} {print}' chr/$id.names - |
  # (samtools faidx для несжатого CS v1.0 пишет в заголовок только первое слово — поэтому ключ и по первому слову)
    $E/bgzip -@ 4 -c > chr/$id.fa.gz && $E/samtools faidx chr/$id.fa.gz && touch chr/$id.done
  echo "$id $(cut -f1 chr/$id.fa.gz.fai | tr '\n' ' ')"; }
export -f one; export E
awk -F'\t' '$1!="CS_v2.1"' assemblies.tsv | xargs -P 4 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch chr/ALL.done
