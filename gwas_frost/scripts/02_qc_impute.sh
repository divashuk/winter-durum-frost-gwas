#!/usr/bin/env bash
# Этап 2. QC генотипов и импутация Beagle 5.4 — отдельно для каждого потока.
# Использование: bash 02_qc_impute.sh stream1|stream2
# Порядок: выбор образцов -> биаллельные SNP -> пропуски по образцу <= 0.3 -> call rate по сайту >= 0.8,
#          MAF >= 0.05, MAC >= 5 -> гетерозиготность по образцу (пометка) -> IBS-дубликаты (>= 0.98, список)
#          -> Beagle 5.4 -> повторный фильтр MAF/MAC -> bed/vcf + числовая матрица для GAPIT.
source "$(dirname "$0")/00_env.sh"
S=$1
VCF="$PROJ/FINAL-181124.vcf.gz"
SAMP="$OUT/tables/02_samples_${S}.txt"
D="$WORK/qc/$S"; mkdir -p "$D"; cd "$D"
LOG="$WORK/logs/02_qc_${S}.log"; : > "$LOG"
say(){ echo "[$(date +%T)] $*" | tee -a "$LOG"; }

say "== $S: образцов в списке $(wc -l < "$SAMP")"
bcftools view -S "$SAMP" --force-samples -m2 -M2 -v snps "$VCF" -Ou \
 | bcftools annotate --rename-chrs "$OUT/tables/02_chr_rename.txt" -Ou \
 | bcftools view -e 'N_PASS(GT="alt")==0 || N_PASS(GT="ref")==0 && N_PASS(GT="het")==0' -Oz -o s0.vcf.gz
bcftools index -f s0.vcf.gz
say "биаллельных полиморфных SNP в подвыборке: $(bcftools index -n s0.vcf.gz); образцов: $(bcftools query -l s0.vcf.gz | wc -l)"

plink2 --vcf s0.vcf.gz --allow-extra-chr --chr-set 14 no-xy --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 50 \
       --make-pgen --out s0 >> "$LOG" 2>&1
# 1) пропуски по образцу
plink2 --pfile s0 --allow-extra-chr --missing --out s0 >> "$LOG" 2>&1
awk 'NR>1 && $NF>0.3{print $1}' s0.smiss > excl_miss.txt
say "образцов с пропусками > 0.3: $(wc -l < excl_miss.txt)"
# 2) call rate сайта, MAF, MAC
plink2 --pfile s0 --allow-extra-chr --remove excl_miss.txt --geno 0.2 --maf 0.05 --mac 5 --make-pgen --out s1 >> "$LOG" 2>&1
say "после geno<=0.2, MAF>=0.05, MAC>=5: $(($(wc -l < s1.pvar)-1)) SNP, $(($(wc -l < s1.psam)-1)) образцов"
plink2 --pfile s1 --allow-extra-chr --missing --freq --out s1 >> "$LOG" 2>&1
# 3) гетерозиготность по образцу (доля гетерозигот среди генотипированных сайтов)
plink2 --pfile s1 --allow-extra-chr --sample-counts --out s1 >> "$LOG" 2>&1
# 4) IBS (plink1.9 --distance ibs, по неимпутированным генотипам)
plink2 --pfile s1 --allow-extra-chr --make-bed --out s1b >> "$LOG" 2>&1
plink --bfile s1b --allow-extra-chr --chr-set 14 --distance square ibs flat-missing --out s1_ibs >> "$LOG" 2>&1
# 5) Beagle 5.4
plink2 --pfile s1 --allow-extra-chr --export vcf bgz id-paste=iid --out s1 >> "$LOG" 2>&1
bcftools index -f s1.vcf.gz
BJAR=$(ls "$ENVBIN"/../share/beagle-5.4*/beagle*.jar | head -1)
say "Beagle: $BJAR"
java -Xmx12g -jar "$BJAR" gt=s1.vcf.gz out=s1_imp nthreads=$THREADS seed=20261025 >> "$LOG" 2>&1
bcftools index -f s1_imp.vcf.gz
# 6) повторный MAF/MAC после импутации
plink2 --vcf s1_imp.vcf.gz --allow-extra-chr --chr-set 14 no-xy --maf 0.05 --mac 5 --make-bed --out final >> "$LOG" 2>&1
plink2 --bfile final --allow-extra-chr --chr-set 14 no-xy --freq --out final >> "$LOG" 2>&1
plink2 --bfile final --allow-extra-chr --chr-set 14 no-xy --export vcf bgz --out final >> "$LOG" 2>&1
plink2 --bfile final --allow-extra-chr --chr-set 14 no-xy --export A --out final >> "$LOG" 2>&1
say "ИТОГ $S: $(wc -l < final.bim) SNP, $(wc -l < final.fam) образцов"
cut -f1 final.bim | sort | uniq -c | awk '{printf "%s:%s ",$2,$1}' | tee -a "$LOG"; echo | tee -a "$LOG"
cp "$LOG" "$OUT/tables/02_qc_${S}.log"
