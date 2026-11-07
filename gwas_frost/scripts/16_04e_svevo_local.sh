#!/usr/bin/env bash
# Раздел 16, шаг 4д (локально, WSL). Svevo v1: хромосомы 5A и 5B уже скачаны в разделе 11 (ref11/svevo_5A.fa, svevo_5B.fa);
# для генов 7B и 1B — только участки ±20 кб вокруг гена (Ensembl REST /sequence/region), целые хромосомы не скачиваются.
#  1) minimap2 -x asm20 запросов (s16/q/queries.fa) -> s16/svevo/svevo.paf (5A, 5B и три участка);
#  2) blastn-short праймеров Pearce -> s16/svevo/pcr.tsv (базы primers11/svevo_5A, svevo_5B раздела 11).
source "$(dirname "$0")/00_env.sh"
O=$WORK/s16/svevo; mkdir -p $O; cd $O
R=https://rest.ensembl.org/sequence/region/triticum_turgidum
reg() { [ -s $1.fa ] || curl -s -m 300 "$R/$2?content-type=text/x-fasta" | sed "1s/.*/>$1/" > $1.fa; }
reg svevo_7B_336600994 7B:336600994..336646124:1      # TRITD7Bv1G107490 ±20 кб
reg svevo_7B_300043888 7B:300043888..300084496:1      # TRITD7Bv1G098240 ±20 кб
reg svevo_1B_625869073 1B:625869073..625912903:1      # TRITD1Bv1G205280 ±20 кб
cat svevo_7B_*.fa svevo_1B_*.fa > svevo_regions.fa; grep -c '>' svevo_regions.fa
minimap2 -t 8 -x asm20 -c --secondary=yes -N 50 -p 0.1 <(cat $WORK/ref11/svevo_5A.fa $WORK/ref11/svevo_5B.fa svevo_regions.fa) $WORK/s16/q/queries.fa > svevo.paf 2> mm.log
: > pcr.tsv
for c in 5A 5B; do blastn -task blastn-short -word_size 7 -evalue 1000 -max_target_seqs 500 -max_hsps 5000 -num_threads 8 \
  -query $WORK/primers11/primers.fa -db $WORK/primers11/svevo_$c \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq" >> pcr.tsv; done
wc -l svevo.paf pcr.tsv
