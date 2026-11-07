#!/usr/bin/env bash
# Раздел 16, шаг 5а. MAFFT (--auto) по каждому гену: genes16/<ген>_{genomic,cds,protein}.fa -> genes16/aln/<ген>_<тип>.aln.fa;
# белки + ортологи ячменя, риса, Brachypodium (Ensembl Compara, s16/compara) -> genes16/aln/<ген>_protein_orth.aln.fa.
# Файлы с одной последовательностью не выравниваются (копируются как есть).
source "$(dirname "$0")/00_env.sh"; set +e   # grep -c = 0 не должен останавливать цикл
MAFFT=/opt/micromamba/envs/s16/bin/mafft
G="$OUT/genes16"; A="$G/aln"; mkdir -p "$A"
cut -f1,2 $WORK/s16/q/candidate_names.tsv | tail -n +2 > $WORK/s16/q/cand_gene.tsv
for f in "$G"/*_genomic.fa "$G"/*_cds.fa "$G"/*_protein.fa; do
  b=$(basename "$f" .fa); n=$(grep -c '>' "$f")
  if [ "$n" -ge 2 ]; then $MAFFT --auto --quiet --thread 6 "$f" > "$A/$b.aln.fa"; elif [ "$n" -eq 1 ]; then cp "$f" "$A/$b.aln.fa"; fi
done
# ортологи: s16/q/file_base.tsv (база файла -> cand_id) пишет 16_04k_integrate.py
while IFS=$'\t' read base cand; do
  f="$G/${base}_protein.fa"; [ -s "$f" ] || continue
  { cat "$f"; seqkit grep -r -p "^${cand//./\\.}\|" $WORK/s16/compara/orthologs_pep.fa | seqkit replace -p '^[^|]+\|' -r ''; } > $WORK/s16/tmp_orth.fa
  $MAFFT --auto --quiet --thread 6 $WORK/s16/tmp_orth.fa > "$A/${base}_protein_orth.aln.fa"
done < $WORK/s16/q/file_base.tsv
ls "$A" | wc -l
