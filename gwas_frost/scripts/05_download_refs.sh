#!/usr/bin/env bash
# Этап 5.0. Загрузка референса Langdon (GCA_051168695.1) и аннотации Svevo v1 (Ensembl Plants) — в WSL, вне OneDrive.
# Аннотации генов для Langdon нет (NCBI, GWH), поэтому гены в интервалах переносятся с Svevo v1 по выравниванию кДНК.
source "$(dirname "$0")/00_env.sh"
R="$WORK/ref"; mkdir -p "$R"; cd "$R"
wget -q -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/168/695/GCA_051168695.1_HAU_Langdon/GCA_051168695.1_HAU_Langdon_genomic.fna.gz
wget -q -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/168/695/GCA_051168695.1_HAU_Langdon/md5checksums.txt -O langdon_md5.txt
grep "_genomic.fna.gz" langdon_md5.txt | grep -v "rna\|cds" | sed 's#\./##' | md5sum -c - | tee md5_check.txt
# только 14 хромосом, переименование CM1182xx.1 -> 1A..7B
zcat GCA_051168695.1_HAU_Langdon_genomic.fna.gz | awk 'BEGIN{split("1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B",c," ")}
  /^>/{keep=0; if(match($0,/^>CM1182([6-7][0-9])\.1/)){i=substr($0,8,2)-59; if(i>=1&&i<=14){keep=1; print ">"c[i]; next}}} keep' > langdon_chr.fa
samtools faidx langdon_chr.fa
cut -f1,2 langdon_chr.fa.fai
wget -q -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/triticum_turgidum/cdna/Triticum_turgidum.Svevo.v1.cdna.all.fa.gz
for c in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B; do
  wget -q -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gff3/triticum_turgidum/Triticum_turgidum.Svevo.v1.63.chromosome.$c.gff3.gz
done
ls -la
