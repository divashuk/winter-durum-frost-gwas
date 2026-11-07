#!/usr/bin/env bash
# Этап 10.2. Гены Svevo v1 в широких интервалах локусов (tables/10_loci.csv) — тот же метод, что на этапе 5:
#  кДНК Svevo той же хромосомы -> minimap2 -x splice:hq на последовательность интервала Langdon;
#  фильтр (идентичность >= 0.98, покрытие транскрипта >= 0.90) и лучший хит на ген — в 10_genes.R.
#  Интервалы обрабатываются по одному (проверка памяти перед каждым).
source "$(dirname "$0")/00_env.sh"
R="$WORK/ref"; G="$WORK/genes10"; mkdir -p "$G"; cd "$G"
"$ENVBIN/Rscript" -e 'q <- data.table::fread(file.path(Sys.getenv("OUT"), "tables/10_loci.csv"));
  writeLines(sprintf("%s:%.0f-%.0f\t%s\t%s", q$chr, q$qtl_start, q$qtl_end, q$locus_id, q$chr), "regions.tsv")'
cat regions.tsv
: > hits.paf
while IFS=$'\t' read -r reg loc chr; do
  echo "$loc $reg $(free -g | awk 'NR==2{print "free_GB="$7}')"
  samtools faidx "$R/langdon_chr.fa" "$reg" | sed "1s/.*/>${loc}|${reg}/" > int.fa
  if [ ! -s cdna_$chr.fa ]; then
    if [ -s "$WORK/genes/cdna_$chr.fa" ]; then cp "$WORK/genes/cdna_$chr.fa" .; else seqkit grep -n -r -p "chromosome:Svevo.v1:${chr}:" "$R/Triticum_turgidum.Svevo.v1.cdna.all.fa.gz" > cdna_$chr.fa; fi
  fi
  minimap2 -t 4 -x splice:hq -c --secondary=no int.fa cdna_$chr.fa 2>/dev/null | cut -f1-12 >> hits.paf
done < regions.tsv
rm -f int.fa
[ -s cdna_headers.txt ] || cp "$WORK/genes/cdna_headers.txt" .
for c in $(cut -f3 regions.tsv | sort -u); do zcat "$R/Triticum_turgidum.Svevo.v1.63.chromosome.$c.gff3.gz" | awk -F'\t' '$3=="gene"'; done > genes.gff3
wc -l hits.paf genes.gff3
