#!/usr/bin/env bash
# Этап 9.6. Транслокация 4AL/5AL/7BS по генам: кДНК всех генов Svevo v1 4A (Ensembl Plants 63) выравниваются
# (minimap2 splice:hq, по одной хромосоме-мишени: 4B, 5B, 7B Langdon); для каждого гена — лучший гомеолог.
# Ожидание: исходный 4AL — гомеологи на 4B; транслоцированный сегмент 5AL — на 5B; 7BS — на 7B.
# Положение генов Q4A.5 на Svevo берётся из tables/05_genes_in_qtl.csv (перенос на Langdon, этап 5).
source "$(dirname "$0")/00_env.sh"
R="$WORK/ref"; T="$WORK/transloc"; mkdir -p "$T"; cd "$T"
seqkit grep -n -r -p "chromosome:Svevo.v1:4A:" "$R/Triticum_turgidum.Svevo.v1.cdna.all.fa.gz" > cdna4A.fa
echo "транскриптов 4A: $(grep -c '>' cdna4A.fa)"
: > genes_hits.paf
for tgt in 4B 5B 7B; do
  free -g | head -2 | tail -1
  samtools faidx "$R/langdon_chr.fa" $tgt > tgt.fa
  minimap2 -t 8 -x splice:hq --secondary=no -c tgt.fa cdna4A.fa 2>> mm_genes.log | cut -f1-12 >> genes_hits.paf
  rm -f tgt.fa
done
grep '>' cdna4A.fa | sed 's/^>//' > cdna4A_headers.txt
"$ENVBIN/Rscript" "$(dirname "$0")/09_4A_gene_synteny.R"
