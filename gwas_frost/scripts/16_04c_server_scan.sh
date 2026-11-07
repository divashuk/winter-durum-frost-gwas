#!/bin/bash
# Раздел 16, шаг 4в (сервер compute-host). По каждой сборке:
#  (память: индекс частями по 4 Гб с --split-prefix, ~7–8 ГБ на процесс; не больше 2 сборок сразу;
#   перед запуском ждёт, пока на сервере доступно ≥ 20 ГБ — на ВМ идут чужие задачи)
#  1) minimap2 -x asm20 запросов (гены-кандидаты Langdon +2 кб/+1 кб и якоря синтении) по целевым хромосомам
#     (Chili, Mahmoudi — по всей сборке): -c --secondary=yes -N 50 -p 0.1 -> scan/<id>.paf;
#  2) BLAST-база хромосом группы 5 (Chili, Mahmoudi — вся сборка) -> blast/<id>;
#  3) blastn-short праймеров Pearce et al. 2013 (Table S1) -> pcr/<id>.tsv (word_size 7, evalue 1000).
# Входы: /work/remote/s16/q/queries.fa, q/primers.fa. Флаги готовности: scan/<id>.done, pcr/<id>.done.
E=/opt/micromamba/envs/st_tools/bin
W=/work/remote/s16; cd $W; mkdir -p scan blast pcr logs
if [ "${1:-}" != run ]; then nohup setsid bash $W/16_04c_server_scan.sh run > logs/scan.log 2>&1 < /dev/null & echo started; exit 0; fi
waitmem() { while [ $(free -g | awk '/Mem/{print $7}') -lt 20 ]; do sleep 60; done; }
one() { id=$1; f=$2; while [ ! -e chr/$id.done ]; do sleep 30; done   # ждёт только своей сборки
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then T="$f"; else T=chr/$id.fa.gz; fi
  if [ ! -e scan/$id.done ]; then
    waitmem; $E/minimap2 -t 8 -x asm20 -I 4G --split-prefix scan/$id.tmp -c --secondary=yes -N 50 -p 0.1 "$T" q/queries.fa > scan/$id.paf 2> logs/mm_$id.log && touch scan/$id.done; fi
  if [ ! -e pcr/$id.done ]; then
    waitmem
    if [ $id = Chili ] || [ $id = Mahmoudi ]; then zcat "$f"; else $E/samtools faidx $T $(cut -f1 $T.fai | grep -E '^5[ABD]$'); fi |
      $E/makeblastdb -in - -dbtype nucl -out blast/$id -title $id -parse_seqids > logs/mkdb_$id.log 2>&1
    $E/blastn -task blastn-short -word_size 7 -evalue 1000 -max_target_seqs 500 -max_hsps 5000 -num_threads 8 -query q/primers.fa -db blast/$id \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq" > pcr/$id.tsv 2> logs/pcr_$id.log && touch pcr/$id.done
  fi; }
export -f one waitmem; export E
awk -F'\t' '$1!="CS_v2.1"' assemblies.tsv | xargs -P 2 -L 1 bash -c 'nice -n 10 bash -c "one $0 $1"'
touch scan/ALL.done
