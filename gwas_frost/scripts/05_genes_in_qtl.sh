#!/usr/bin/env bash
# Этап 5.3. Гены в интервалах QTL. У сборки Langdon нет аннотации генов, поэтому:
#  1) вырезаем последовательности интервалов Langdon (tables/05_qtl_intervals.csv);
#  2) выравниваем на них кДНК Svevo v1 (Ensembl Plants 63) той же хромосомы: minimap2 -x splice:hq;
#  3) оставляем выравнивания с идентичностью >= 0.98 и покрытием транскрипта >= 0.90, лучший хит на ген;
#  4) описание гена — только из аннотации Svevo (поле description: в заголовке кДНК / GFF3), без собственных догадок.
source "$(dirname "$0")/00_env.sh"
R="$WORK/ref"; G="$WORK/genes"; mkdir -p "$G"; cd "$G"
# только воспроизводимые QTL (столбец reliability); единичные сигналы одной модели гены не получают
Rscript -e 'q <- data.table::fread(file.path(Sys.getenv("OUT"), "tables/05_qtl_intervals.csv")); q <- q[reliability == "воспроизводимый"];
  writeLines(sprintf("%s:%.0f-%.0f\t%s\t%s", q$chr, q$start, q$end, q$QTL, q$chr), "regions.tsv")'
cat regions.tsv
: > hits.paf
while IFS=$'\t' read -r reg qtl chr; do
  samtools faidx "$R/langdon_chr.fa" "$reg" | sed "1s/.*/>${qtl}|${reg}/" > int.fa
  [ -s cdna_$chr.fa ] || seqkit grep -n -r -p "chromosome:Svevo.v1:${chr}:" "$R/Triticum_turgidum.Svevo.v1.cdna.all.fa.gz" > cdna_$chr.fa
  minimap2 -t 4 -x splice:hq -c --secondary=no int.fa cdna_$chr.fa 2>/dev/null >> hits.paf
done < regions.tsv
zcat "$R/Triticum_turgidum.Svevo.v1.cdna.all.fa.gz" | grep ">" | sed 's/^>//' > cdna_headers.txt
for c in $(cut -f3 regions.tsv | sort -u); do zcat "$R/Triticum_turgidum.Svevo.v1.63.chromosome.$c.gff3.gz" | awk -F'\t' '$3=="gene"'; done > genes.gff3
Rscript "$(dirname "$0")/05_genes_in_qtl.R"
