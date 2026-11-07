#!/bin/bash
# Раздел 16, шаг 4н (сервер compute-host). BLAST-базы и blastn-short праймеров для сборок, где шаг 3 скрипта 16_04c не выполнился
# (в первом запуске 16_04c ожидание памяти стояло внутри конвейера перед makeblastdb, базы получились пустыми).
#  праймеры Pearce (q/primers.fa) -> pcr/<id>.tsv;  фланки SNP GWAS ±150 п.н. (q/gwas_flanks.fa, blastn) -> pcr/<id>_flanks.tsv
E=/opt/micromamba/envs/st_tools/bin
W=/work/remote/s16; cd $W; mkdir -p blast pcr logs
if [ "${1:-}" != run ]; then nohup setsid bash $W/16_04m_server_pcr.sh run > logs/pcr_all.log 2>&1 < /dev/null & echo started; exit 0; fi
waitmem() { while [ $(free -g | awk '/Mem/{print $7}') -lt 15 ]; do sleep 60; done; }
while IFS=$'\t' read id f; do
  [ "$id" = CS_v2.1 ] && continue; [ -e pcr/$id.done ] && continue
  while [ ! -e chr/$id.done ]; do sleep 60; done
  if [ $id = Chili ] || [ $id = Mahmoudi ]; then T="$f"; else T=chr/$id.fa.gz; fi
  waitmem
  if [ ! -s blast/$id.ndb ] && [ ! -s blast/$id.nal ]; then
    if [ $id = Chili ] || [ $id = Mahmoudi ]; then zcat "$f"; else $E/samtools faidx $T $(cut -f1 $T.fai | grep -E '^5[ABD]$'); fi > blast/$id.tmp.fa
    nice -n 10 $E/makeblastdb -in blast/$id.tmp.fa -dbtype nucl -out blast/$id -title $id -parse_seqids > logs/mkdb_$id.log 2>&1; rm -f blast/$id.tmp.fa
  fi
  nice -n 10 $E/blastn -task blastn-short -word_size 7 -evalue 1000 -max_target_seqs 500 -max_hsps 5000 -num_threads 8 -query q/primers.fa -db blast/$id \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq" > pcr/$id.tsv 2> logs/pcr_$id.log &&
  nice -n 10 $E/blastn -task blastn -evalue 1e-10 -max_target_seqs 50 -max_hsps 20 -num_threads 8 -query q/gwas_flanks.fa -db blast/$id \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq" > pcr/${id}_flanks.tsv 2>> logs/pcr_$id.log &&
  touch pcr/$id.done
done < assemblies.tsv
touch pcr/ALL.done
