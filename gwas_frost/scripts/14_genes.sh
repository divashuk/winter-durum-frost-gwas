#!/usr/bin/env bash
# Этап 14.4 (часть 2). Гены Svevo v1 (Ensembl Plants 63, все гены — HC: logic_name=gff3_hc_genes) в интервалах tables/14_gene_intervals.csv:
#  кДНК генов той же хромосомы Svevo -> minimap2 -x splice:hq на последовательность интервала Langdon (как в разделах 5, 10, 11);
#  фильтр идентичности/покрытия и коллинеарности — в 14_genes.R. Интервалы — по одному, с проверкой памяти.
source "$(dirname "$0")/00_env.sh"
R="$WORK/ref"; G="$WORK/genes14"; mkdir -p "$G"; cd "$G"
"$ENVBIN/Rscript" -e 'q <- data.table::fread(file.path(Sys.getenv("OUT"), "tables/14_gene_intervals.csv"));
  writeLines(sprintf("%s:%.0f-%.0f\t%s\t%s", q$chr, q$start, q$end, q$interval_id, q$chr), "regions.tsv")'
cat regions.tsv
: > hits.paf
while IFS=$'\t' read -r reg iid chr; do
  echo "$iid $reg $(free -g | awk 'NR==2{print "free_GB="$7}')"
  samtools faidx "$R/langdon_chr.fa" "$reg" | sed "1s/.*/>${iid}|${reg}/" > int.fa
  if [ ! -s cdna_$chr.fa ]; then
    if [ -s "$WORK/genes10/cdna_$chr.fa" ]; then cp "$WORK/genes10/cdna_$chr.fa" .; else seqkit grep -n -r -p "chromosome:Svevo.v1:${chr}:" "$R/Triticum_turgidum.Svevo.v1.cdna.all.fa.gz" > cdna_$chr.fa; fi
  fi
  minimap2 -t 8 -x splice:hq -c --secondary=no int.fa cdna_$chr.fa 2>/dev/null | cut -f1-12 >> hits.paf
done < regions.tsv
rm -f int.fa
[ -s cdna_headers.txt ] || cp "$WORK/genes10/cdna_headers.txt" .
for c in $(cut -f3 regions.tsv | sort -u); do zcat "$R/Triticum_turgidum.Svevo.v1.63.chromosome.$c.gff3.gz" | awk -F'\t' '$3=="gene"'; done > genes.gff3
wc -l hits.paf genes.gff3
