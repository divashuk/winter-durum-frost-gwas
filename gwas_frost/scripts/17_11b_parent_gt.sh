#!/usr/bin/env bash
# Раздел 17, шаг 7б (повтор 16_07a). Исходный VCF коллекции FINAL-181124.vcf.gz (только GT; читается потоком, не меняется):
#  1) генотипы родителей F2 и Svevo коллекции (TP_MO_216) на SNP типа (б) (s17/kasp/snps_b.tsv) и на всех SNP внутри запросов генов
#     части A (ген + 2 кб + 1 кб, s17/q/candidates.tsv) -> s17/kasp/parents_gt.tsv;
#  2) все SNP обоих исходных VCF (коллекция и F2) в ± 150 п.н. от SNP типа (б) -> s17/kasp/gwas_flank_vcf_snps.tsv (маска под праймеры);
#  3) шаблоны Langdon ± 150 п.н. -> s17/kasp/gwas_flanks.fa (имя SNP_<хромосома>_<позиция>);
#  4) проверка на паралоги (ведущие SNP 1B, 7B, 2B): доля гетерозигот по SNP у образцов потока 2 (tables/02_samples_stream2.txt,
#     326 образцов коллекции) и средняя доля гетерозигот по всем SNP VCF у тех же образцов -> s17/kasp/het_by_snp.tsv, het_genome.tsv.
source "$(dirname "$0")/00_env.sh"
K=$WORK/s17/kasp; mkdir -p $K; cd $K
S="TP_MO_081,TP_MO_044,TP_MO_054,TP_MO_014,TP_MO_005,TP_MO_037,TP_MO_171,TP_MO_190,TP_MO_196,TP_MO_302,TP_MO_304,TP_MO_350,TP_MO_216"
awk -F'\t' 'NR==FNR{split($0,a," "); m[a[2]]=a[1]; next} FNR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
  {print m[$h["chr"]] "\t" $h["qstart"] "\t" $h["qend"]}' "$OUT/tables/02_chr_rename.txt" $WORK/s17/q/candidates.tsv > reg.tsv
awk -F'\t' 'NR==FNR{split($0,a," "); m[a[2]]=a[1]; next} FNR==1{next} {print m[$4] "\t" $5 "\t" $5}' "$OUT/tables/02_chr_rename.txt" snps_b.tsv >> reg.tsv
sort -k1,1 -k2,2n reg.tsv > reg.s.tsv
bcftools view -s $S -T reg.s.tsv "$PROJ/FINAL-181124.vcf.gz" -Ou | bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > parents_gt.tsv
wc -l parents_gt.tsv
awk -F'\t' 'NR==FNR{split($0,a," "); m[a[2]]=a[1]; next} FNR==1{next} {print m[$4] "\t" $5-150 "\t" $5+150}' "$OUT/tables/02_chr_rename.txt" snps_b.tsv | sort -k1,1 -k2,2n > flank_reg.tsv
for v in FINAL-181124 FINAL-251124; do bcftools view -G -T flank_reg.tsv "$PROJ/$v.vcf.gz" -Ou | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t$v\n"; done > gwas_flank_vcf_snps.tsv
cat gwas_flank_vcf_snps.tsv
: > gwas_flanks.fa
tail -n +2 snps_b.tsv | while IFS=$'\t' read snp locus role c p; do
  samtools faidx $WORK/primers11/langdon_$c.fa langdon_$c:$((p-150))-$((p+150)) | sed "1s/.*/>SNP_${c}_${p}/" >> gwas_flanks.fa; done
grep -c '>' gwas_flanks.fa
# доля гетерозигот: образцы потока 2
tr -d '\r' < "$OUT/tables/02_samples_stream2.txt" > s2.txt
bcftools view -S s2.txt --force-samples "$PROJ/FINAL-181124.vcf.gz" -Ou | bcftools query -f '%CHROM\t%POS[\t%GT]\n' |
  awk -F'\t' 'NR==FNR{split($0,a," "); m[a[1]]=a[2]; next}
    { h=0; n=0; for(i=3;i<=NF;i++){ if($i!="./." && $i!="."){ n++; split($i,g,/[\/|]/); if(g[1]!=g[2]) h++ } } c=m[$1]; print c "\t" $2 "\t" n "\t" h "\t" (n? h/n : "NA") }' \
    "$OUT/tables/02_chr_rename.txt" - > het_all.tsv
awk -F'\t' '$5!="NA"{s+=$5; k++} END{printf "mean_het_all_snps\t%.5f\t%d\n", s/k, k}' het_all.tsv > het_genome.tsv
awk -F'\t' 'NR==FNR{if(FNR>1) w[$4":"$5]=1; next} ($1":"$2 in w)' snps_b.tsv het_all.tsv > het_by_snp.tsv
cat het_genome.tsv het_by_snp.tsv
