#!/bin/bash
# Раздел 17, шаг 4е (сервер compute-host; повтор 16_04i + 16_04l). Второй проход по списку участков q/pass2_<id>.tsv (17_08_pass2_regions.py):
#  win  — синтенные окна -> p2/<id>_win.fa; miniprot белков генов части A (q/prot_long.fa: Svevo + ортолог CS) с параметрами
#         окончательного прогона раздела 16 (16_04l_local_pass2.sh): --gff --trans -G 20k -p 0.3 --outs=0.3 -N 20 -> p2/<id>_win_long.gff;
#         белков CBF нет — p2/<id>_win.gff пустой;
#  copy — ген + 2 кб + 1 кб каждой копии в ориентации гена -> p2/<id>_copies.fa.
# Хромосомы: группы 1, 5, 7 — /work/remote/s16/chr/<id>.fa.gz; 2, 4 — chr/<id>_g24.fa.gz; Chili, Mahmoudi — исходный файл (индекс раздела 16).
E=/opt/micromamba/envs/st_tools/bin; MP=/opt/micromamba/envs/st_age/bin/miniprot
W=/work/remote/s17; W16=/work/remote/s16; cd $W; mkdir -p p2 logs
if [ "${1:-}" != run ]; then nohup setsid bash $W/17_08b_server_pass2.sh run > logs/pass2.log 2>&1 < /dev/null & echo started; exit 0; fi
fx() { id=$1; f=$2; reg=$3; shift 3      # faidx участка из файла, где есть его хромосома
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then $E/samtools faidx --fai-idx $W16/chr/$id.fai --gzi-idx $W16/chr/$id.gzi "$f" "$@" "$reg"; return; fi
  c=${reg%:*}; for t in $W16/chr/$id.fa.gz chr/${id}_g24.fa.gz; do
    if cut -f1 $t.fai | grep -qx "$c"; then $E/samtools faidx $t "$@" "$reg"; return; fi; done; }
one() { id=$1; f=$2; L=q/pass2_$id.tsv; [ -s $L ] || return 0; [ -e p2/$id.done ] && return 0
  : > p2/${id}_win.fa
  awk -F'\t' '$1=="win"{print $3":"$4"-"$5}' $L | while read r; do fx $id "$f" $r >> p2/${id}_win.fa; done
  : > p2/${id}_win.gff
  $MP -t 6 --gff --trans -G 20k -p 0.3 --outs=0.3 -N 20 p2/${id}_win.fa q/prot_long.fa > p2/${id}_win_long.gff 2> logs/mp_$id.log
  : > p2/${id}_copies.fa
  awk -F'\t' '$1=="copy"' $L | while IFS=$'\t' read k name c s e st; do
    if [ "$st" = "-" ]; then fx $id "$f" "$c:$s-$e" -i; else fx $id "$f" "$c:$s-$e"; fi | sed "1s/.*/>$name/" >> p2/${id}_copies.fa; done
  touch p2/$id.done; echo "$id $(date)"; }
export -f one fx; export E MP W W16
awk -F'\t' '$1!="CS_v2.1" && $1!="Svevo"' $W16/assemblies.tsv | xargs -P 4 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch p2/ALL.done
