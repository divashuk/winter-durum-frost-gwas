#!/usr/bin/env bash
# Раздел 17, часть A, шаг 2. Перенос генов Svevo v1 на новые интервалы I09–I12 (tables/17_gene_intervals.csv) и аннотация —
#  скрипты раздела 14 без изменения метода: 14_genes.sh (minimap2 splice:hq кДНК той же хромосомы на интервал Langdon),
#  14_genes.R (идентичность ≥ 0.98, покрытие ≥ 0.90, коллинеарность ± 5 Мб от скользящей медианы, r² ближайшего SNP),
#  14_biomart.py (InterPro, GO, Compara). Меняются только пути: genes14 -> genes17, список интервалов -> genes17/intervals_new.csv,
#  счётчик генов -> genes17/gene_counts.csv. Кэш BioMart раздела 14 копируется (дозапрашиваются новые гены).
source "$(dirname "$0")/00_env.sh"
S="$OUT/scripts"; G=$WORK/genes17; mkdir -p $G
"$ENVBIN/Rscript" -e 'x <- data.table::fread(file.path(Sys.getenv("OUT"), "tables/17_gene_intervals.csv")); data.table::fwrite(x[interval_id %in% c("I09","I10","I11","I12")], "/work/scratch/genes17/intervals_new.csv")'
for f in $WORK/genes14/biomart_*.tsv; do [ -e $G/$(basename $f) ] || cp $f $G/; done
sed -e 's#genes14#genes17#g' -e 's#file.path(Sys.getenv("OUT"), "tables/14_gene_intervals.csv")#"/work/scratch/genes17/intervals_new.csv"#' \
    "$S/14_genes.sh" > $G/run_14_genes.sh
sed -e 's#"genes14"#"genes17"#' -e 's#file.path(OUT, "tables/14_gene_intervals.csv")#"/work/scratch/genes17/intervals_new.csv"#' \
    -e 's#file.path(OUT, "tables/14_gene_counts_by_interval.csv")#file.path(G, "gene_counts.csv")#' "$S/14_genes.R" > $G/run_14_genes.R
cp "$S/00_env.sh" $G/00_env.sh
bash $G/run_14_genes.sh
"$ENVBIN/Rscript" $G/run_14_genes.R
