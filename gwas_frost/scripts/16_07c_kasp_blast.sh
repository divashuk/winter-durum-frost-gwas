#!/usr/bin/env bash
# Раздел 16, шаг 7в. BLAST для KASP:
#   режим tmpl   — шаблоны (Langdon ± 150 п.н.; s16/kasp/templates.fa), blastn, evalue 1e-10:
#                  по 14 хромосомам Langdon (базы primers11/langdon_*), Svevo 5A/5B (primers11/svevo_*) и CS RefSeq v1.0 целиком (сервер);
#                  нецелевые попадания = гомеологи и паралоги, по ним выбираются позиции, отличающие целевую копию;
#   режим primer — праймеры (специфическая часть без хвоста; s16/kasp/primers_cand.fa), blastn-short (word_size 7, evalue 1000)
#                  по тем же базам. Svevo — только 5A и 5B (остальные хромосомы Svevo не скачивались).
# Выход: s16/kasp/<режим>_{langdon,svevo,cs}.tsv
source "$(dirname "$0")/00_env.sh"; set +e
K=$WORK/s16/kasp; P=$WORK/primers11; cd $K
mode=$1
if [ "$mode" = tmpl ]; then Q=templates.fa; OPT="-task blastn -evalue 1e-10 -max_target_seqs 50 -max_hsps 50"; FLT='1'
else Q=primers_cand.fa; OPT="-task blastn-short -word_size 7 -evalue 1000 -max_target_seqs 500 -max_hsps 2000"
  FLT='$8==$13 && $6==0 && ($5+$7-1)<=2'; fi     # праймеры: сразу только попадания до 3′-конца, без гэпов, ≤ 2 несовпадений
F="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq"
: "${PART:=all}"     # PART=local — только Langdon и Svevo, PART=cs — только CS (для параллельного запуска)
if [ "$PART" != cs ]; then
: > ${mode}_langdon.tsv; : > ${mode}_svevo.tsv
for c in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B; do blastn $OPT -num_threads 8 -query $Q -db $P/langdon_$c -outfmt "$F" | awk -F'	' "$FLT" >> ${mode}_langdon.tsv; done
for c in 5A 5B; do blastn $OPT -num_threads 8 -query $Q -db $P/svevo_$c -outfmt "$F" | awk -F'	' "$FLT" >> ${mode}_svevo.tsv; done
fi
[ "$PART" != local ] && ssh -o BatchMode=yes compute-host "cat > /tmp/kasp_q_$mode.fa; nice -n 10 /opt/micromamba/envs/st_tools/bin/blastn $OPT -num_threads 8 -query /tmp/kasp_q_$mode.fa -db /work/remote/reference/CS_RefSeq1.0 -outfmt '$F' | awk -F'\t' '$FLT'" < $Q > ${mode}_cs.tsv
wc -l ${mode}_*.tsv
