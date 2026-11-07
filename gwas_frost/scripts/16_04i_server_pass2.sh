#!/bin/bash
# Раздел 16, шаг 4и (сервер compute-host). Второй проход по каждой сборке (список участков q/pass2_<id>.tsv из 16_04g_pass2_regions.py):
#  win  — синтенные окна групп (якоря ± 1 Мб) -> p2/<id>_win.fa (samtools faidx);
#         miniprot --gff --trans -G 3k (белки кандидатов: CBF CS IPR045277 гр. 5, T. monococcum, Svevo; белки прочих кандидатов Svevo/CS)
#         -> p2/<id>_win.gff  (интрон до 3 кб достаточен для CBF; для многоэкзонных генов — второй запуск с -G 20k -> p2/<id>_win_long.gff);
#  copy — последовательность каждой копии (ген + 2 кб + 1 кб в ориентации гена) -> p2/<id>_copies.fa.
# Chili, Mahmoudi: faidx по исходному файлу хранилища с индексом chr/<id>.fai/.gzi.
E=/opt/micromamba/envs/st_tools/bin; MP=/opt/micromamba/envs/st_age/bin/miniprot
W=/work/remote/s16; cd $W; mkdir -p p2 logs
if [ "${1:-}" != run ]; then nohup setsid bash $W/16_04i_server_pass2.sh run > logs/pass2.log 2>&1 < /dev/null & echo started; exit 0; fi
one() { id=$1; f=$2; L=q/pass2_$id.tsv; [ -s $L ] || return 0; [ -e p2/$id.done ] && return 0
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then FX="$E/samtools faidx --fai-idx chr/$id.fai --gzi-idx chr/$id.gzi $f"; else FX="$E/samtools faidx chr/$id.fa.gz"; fi
  awk -F'\t' '$1=="win"{print $3":"$4"-"$5}' $L > p2/$id.win.reg
  $FX -r p2/$id.win.reg > p2/${id}_win.fa
  $MP -t 6 --gff --trans -G 3k -p 0.3 --outs=0.3 -N 100 --outn=1000 p2/${id}_win.fa q/prot2.fa > p2/${id}_win.gff 2> logs/mp_$id.log
  $MP -t 6 --gff --trans -G 20k -N 5 p2/${id}_win.fa q/prot_long.fa > p2/${id}_win_long.gff 2>> logs/mp_$id.log
  : > p2/${id}_copies.fa
  awk -F'\t' '$1=="copy"' $L | while IFS=$'\t' read k name c s e st; do
    if [ "$st" = "-" ]; then $FX -i "$c:$s-$e"; else $FX "$c:$s-$e"; fi | sed "1s/.*/>$name/" >> p2/${id}_copies.fa; done
  touch p2/$id.done; }
export -f one; export E MP
awk -F'\t' '$1!="CS_v2.1"' assemblies.tsv | xargs -P 4 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch p2/ALL.done
