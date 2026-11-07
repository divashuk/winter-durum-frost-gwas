#!/usr/bin/env bash
# Раздел 17, шаг 5а. MAFFT --auto по генам части A — scripts/16_05a_align.sh без изменения метода (замены: genes16 -> genes17, s16 -> s17).
# Для рисунка выравниваний (17_13c) в genes17/aln копируются выравнивания TRITD1Bv1G205280 раздела 16.
source "$(dirname "$0")/00_env.sh"
sed -e 's#G="$OUT/genes16"#G="$OUT/genes17"#' -e 's#\$WORK/s16/#$WORK/s17/#g' "$OUT/scripts/16_05a_align.sh" > $WORK/s17/run_16_05a.sh
grep -c 's17' $WORK/s17/run_16_05a.sh
cp "$OUT/scripts/00_env.sh" $WORK/s17/00_env.sh
bash $WORK/s17/run_16_05a.sh
cp "$OUT"/genes16/aln/TRITD1Bv1G205280_* "$OUT/genes17/aln/"
