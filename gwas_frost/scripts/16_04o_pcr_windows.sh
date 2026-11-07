#!/usr/bin/env bash
# Раздел 16, шаг 4п (локально). In silico ПЦР праймеров Pearce по синтенным окнам каждой сборки (p2/<g>_win.fa: окна групп
# на любых хромосомах, в т.ч. 7B у сборок с транслокацией, скаффолды Chili/Mahmoudi). blastn-short как в 16_04c -> s16/pcr_win/<g>.tsv
source "$(dirname "$0")/00_env.sh"; set +e
S=$WORK/s16; O=$S/pcr_win; mkdir -p $O; cd $O
for w in $S/p2/*_win.fa; do g=$(basename $w _win.fa)
  makeblastdb -in $w -dbtype nucl -out db_$g > /dev/null
  blastn -task blastn-short -word_size 7 -evalue 1000 -max_target_seqs 500 -max_hsps 5000 -num_threads 8 -query $WORK/primers11/primers.fa -db db_$g \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq" |
    awk -F'\t' '$8==$13 && $6==0 && ($5+$7-1)<=2' > $g.tsv
  echo "$g $(wc -l < $g.tsv)"
done
