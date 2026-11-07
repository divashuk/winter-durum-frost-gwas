#!/usr/bin/env bash
# Раздел 17, шаг 7г. BLAST шаблонов (режим tmpl) и праймеров (режим primer) — scripts/16_07c_kasp_blast.sh без изменения метода:
#  Langdon — 14 хромосом (primers11/langdon_*), CS RefSeq v1.0 целиком (сервер). Svevo: в 16 — только 5A и 5B (primers11/svevo_*);
#  здесь ещё база участков Svevo раздела 17 (s17/svevo/svevo_regions.fa: синтенные участки генов части A, 3–6 Мб каждый).
# Запуск: bash scripts/17_14b_kasp_blast.sh tmpl|primer  (PART=local / cs — как в 16)
source "$(dirname "$0")/00_env.sh"
SV=$WORK/s17/svevo; [ -s $SV/db.nsq ] || [ -s $SV/db.00.nsq ] || makeblastdb -in $SV/svevo_regions.fa -dbtype nucl -out $SV/db > /dev/null
sed -e 's#K=$WORK/s16/kasp#K=$WORK/s17/kasp#' \
    -e 's#for c in 5A 5B; do blastn $OPT -num_threads 8 -query $Q -db $P/svevo_$c#for c in $P/svevo_5A $P/svevo_5B $WORK/s17/svevo/db; do blastn $OPT -num_threads 8 -query $Q -db $c#' \
    -e 's#/tmp/kasp_q_$mode.fa#/tmp/kasp17_q_$mode.fa#g' \
    "$OUT/scripts/16_07c_kasp_blast.sh" > $WORK/s17/run_16_07c.sh
grep -c 's17' $WORK/s17/run_16_07c.sh
bash $WORK/s17/run_16_07c.sh "$@"
