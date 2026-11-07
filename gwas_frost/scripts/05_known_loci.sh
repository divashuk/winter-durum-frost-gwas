#!/usr/bin/env bash
# Этап 5.1. Позиции известных локусов и KASP-маркеров на Langdon.
#  - гены: последовательности NCBI (accession в known_genes.tsv) -> minimap2 (-x splice для генов, asm20 для BAC) на хромосомы Langdon
#  - KASP: праймеры (аллель-специфичный без хвоста FAM/HEX + общий обратный) -> blastn-short по ожидаемой хромосоме;
#    ампликон = пара попаданий на разных цепях в пределах 1 кб.
source "$(dirname "$0")/00_env.sh"
R="$WORK/ref"; K="$WORK/known"; mkdir -p "$K"; cd "$K"
cat > known_genes.tsv <<'EOF'
VRN-A1	AY747598	5A	T. turgidum Langdon VRN-A1, полный ген
VRN-B1	AY747602	5B	T. turgidum Langdon VRN-B1, полный ген
VRN-B3_FT-B1	DQ890164	7B	T. turgidum Langdon VRN3 (vrn-B3)
Ppd-A1	AB693038	2A	T. durum Ppd-A1, полный ген
Ppd-B1	DQ885757	2B	T. aestivum CS 2B PRR (Ppd-B1)
Rht-A1(→Rht-B1 гомеолог)	JF930277	4A,4B	T. aestivum Rht-A1a; попадание на 4B = гомеолог Rht-B1
Fr-Am2_BAC	AY951944	5A,5B	T. monococcum BAC 21C6, локус Fr-Am2 (кластер CBF)
CBF12	EU076381	5A,5B	T. monococcum CBF12
CBF15	EU076383	5A,5B	T. monococcum CBF15
EOF
: > known_genes.fa
while IFS=$'\t' read -r name acc chr note; do
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text" \
    | sed "1s/^>.*/>${name}|${acc}/" >> known_genes.fa
  sleep 0.5
done < known_genes.tsv
grep -c ">" known_genes.fa
# гены: minimap2 asm20 ПО ОДНОЙ ожидаемой хромосоме (индекс всего генома 10 Гб требует > 21 Гб ОЗУ — был OOM)
: > known_genes.paf
for c in 1B 2A 2B 4A 4B 5A 5B 7B; do
  [ -s chr_$c.fa ] || samtools faidx "$R/langdon_chr.fa" $c > chr_$c.fa
  awk -F'	' -v c=$c 'index(","$3",", ","c",")>0{print $1"|"$2}' known_genes.tsv > q_$c.txt
  [ -s q_$c.txt ] || continue
  seqkit grep -n -f q_$c.txt known_genes.fa > q_$c.fa
  minimap2 -t 4 -x asm20 -c --secondary=yes -N 5 chr_$c.fa q_$c.fa >> known_genes.paf 2>> mm.log
done
awk 'BEGIN{OFS="\t"; print "query","qlen","qstart","qend","strand","chr","tstart","tend","nmatch","alnlen","mapq","identity","qcov"}
     {print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,sprintf("%.4f",$10/$11),sprintf("%.3f",($4-$3)/$2)}' known_genes.paf > known_genes_hits.tsv
# KASP-праймеры
cat > kasp_primers.fa <<'EOF'
>S2269949_ASfam
ACTGCAGCGAGGATTCCCA
>S2269949_AShex
CTACTGCAGCGAGGATTCCCG
>S2269949_R
AGCCACCTCTGAGATGCCGCTT
>S1077313_ASfam
ATATGTCTCCCTCACCATGGCG
>S1077313_R
GTCAACTGCTTTGTACTGATTCGTGATAT
>S1862541_ASfam
CAGCTCACCTCGCCTCTGAG
>S1862541_R
TGAATTCATGGTCAGCTCGGTCGAA
>S1298957_ASfam
CAGCTGGATATTTCCCCGCTC
>S1298957_R
CCACATCTTGTGGGTTAATTGCGGTT
>S1051014_ASfam
AGCATTTCAGGATGCAAAGGCATC
>S1051014_R
GTAGAAGCAAATTTCGTGGTGAGCGTT
>Mor587_1B_41099587_A
CCGTCCATGTTCTCTCGCAT
>Mor587_1B_41099587_C
CTGGATCCATCACGCCAAG
>Mor534_5B_517276534_A
GGATGAGGATGAGCAAAGCAA
>Mor534_5B_517276534_C
CCATGCATGGGCCACTTTTTC
>Mor866_7B_598228866_A
GGTCCGAATGTTATAAGGCCC
>Mor866_7B_598228866_C
CGAAATCGACCTACAAACTGTTG
EOF
cat chr_5A.fa chr_5B.fa chr_1B.fa chr_7B.fa > sub.fa
makeblastdb -in sub.fa -dbtype nucl -out sub > /dev/null
blastn -task blastn-short -query kasp_primers.fa -db sub -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
       -evalue 1 -max_target_seqs 50 -max_hsps 20 -num_threads 4 -word_size 11 > kasp_blast.tsv
# полные попадания (длина = длине праймера, <=1 несовпадение, 3'-конец включён)
awk '$4>=$13-1 && $5<=1 && $8==$13' kasp_blast.tsv > kasp_full_hits.tsv
wc -l kasp_full_hits.tsv
cp known_genes_hits.tsv "$OUT/tables/05_known_genes_hits_raw.tsv"
cp kasp_full_hits.tsv "$OUT/tables/05_kasp_primer_hits_raw.tsv"
cp known_genes.tsv "$OUT/tables/05_known_genes_queries.tsv"
