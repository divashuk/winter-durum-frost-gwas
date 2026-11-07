#!/usr/bin/env bash
# Этап 10.7. KASP: фланки ±100 п.н. (Langdon), соседние варианты обоих VCF, уникальность выравнивания фланка по геному.
#  Уникальность: 201-п.н. фланк (референсный аллель) -> minimap2 -x sr --secondary=yes -N 50 по каждой хромосоме Langdon отдельно
#  (индекс всего генома не помещается в ОЗУ). Хит засчитывается при длине выравнивания >= 150 п.н.; идентичность считается в 10_kasp_final.R.
source "$(dirname "$0")/00_env.sh"
K="$WORK/kasp"; R="$WORK/ref"; mkdir -p "$K"; cd "$K"
"$ENVBIN/Rscript" "$(dirname "$0")/10_kasp_select.R"
samtools faidx "$R/langdon_chr.fa" -r regions.txt > flanks.fa
# все варианты исходных VCF (до КК); контиги CM118260.1–CM118273.1 = 1A…7B
CH=(1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B)
: > all_variants.tsv
for v in "$PROJ/FINAL-181124.vcf.gz" "$PROJ/FINAL-251124.vcf.gz"; do bcftools query -f '%CHROM\t%POS\n' "$v" >> all_variants.tsv; done
awk -v OFS='\t' 'BEGIN{split("1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B",c," "); for(i=0;i<14;i++) m["CM1182" (60+i) ".1"]=c[i+1]} {print ($1 in m? m[$1]:$1), $2}' all_variants.tsv | sort -u > all_variants_named.tsv
: > flank_hits.paf
for c in "${CH[@]}"; do
  echo "$c $(free -g | awk 'NR==2{print "free_GB="$7}')"
  samtools faidx "$R/langdon_chr.fa" "$c" > chr.fa
  minimap2 -t 8 -x sr --secondary=yes -N 50 -p 0.5 -c chr.fa flanks.fa 2>/dev/null | cut -f1-12 >> flank_hits.paf
  rm -f chr.fa
done
wc -l flank_hits.paf all_variants_named.tsv
"$ENVBIN/Rscript" "$(dirname "$0")/10_kasp_final.R"
