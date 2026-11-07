#!/usr/bin/env bash
# Раздел 16, шаг 7а. Генотипы GBS (исходный VCF коллекции FINAL-181124.vcf.gz, только GT; файл читается потоком, не меняется)
# родителей F2 и Svevo (TP_MO_216) на 4 SNP GWAS и на всех SNP внутри запросов генов-кандидатов (ген + 2 кб + 1 кб).
# Выход: WSL s16/kasp/parents_gt.tsv
source "$(dirname "$0")/00_env.sh"
O=$WORK/s16/kasp; mkdir -p $O; cd $O
S="TP_MO_081,TP_MO_044,TP_MO_054,TP_MO_014,TP_MO_005,TP_MO_037,TP_MO_171,TP_MO_190,TP_MO_196,TP_MO_302,TP_MO_304,TP_MO_350,TP_MO_216"
# регионы в именах контигов VCF (CM1182xx.1); соответствие — tables/02_chr_rename.txt (контиг -> хромосома)
awk -F'\t' 'NR==FNR{split($0,a," "); m[a[2]]=a[1]; next} FNR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
  {print m[$h["chr"]] "\t" $h["qstart"] "\t" $h["qend"]}' "$OUT/tables/02_chr_rename.txt" $WORK/s16/q/candidates.tsv > reg.tsv
awk 'NR==FNR{m[$2]=$1; next} {print m[$1] "\t" $2 "\t" $2}' "$OUT/tables/02_chr_rename.txt" - >> reg.tsv <<'R'
5B	535723924
5B	535618213
5B	522742783
5B	532482598
R
sort -k1,1 -k2,2n reg.tsv > reg.s.tsv
bcftools view -s $S -T reg.s.tsv "$PROJ/FINAL-181124.vcf.gz" -Ou | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > parents_gt.tsv
wc -l parents_gt.tsv; head -3 parents_gt.tsv
# все SNP VCF коллекции (и F2) в ± 150 п.н. от 4 SNP GWAS — позиции под праймеры не ставятся (s16/kasp/gwas_flank_vcf_snps.tsv)
for p in 535723924 535618213 522742783 532482598; do printf 'CM118269.1\t%s\t%s\n' $((p-150)) $((p+150)); done > flank_reg.tsv
for v in FINAL-181124 FINAL-251124; do bcftools view -G -T flank_reg.tsv "$PROJ/$v.vcf.gz" -Ou | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t$v\n"; done > gwas_flank_vcf_snps.tsv
cat gwas_flank_vcf_snps.tsv
