#!/usr/bin/env bash
# Раздел 16, шаг 4м (локально, WSL). Уточнение второго прохода по синтенным окнам (p2/<g>_win.fa, перенесены с сервера — только окна):
#  1) minimap2 -x asm20 дополнительных запросов (s16/q/queries_extra.fa: вторая копия TRITD5Bv1G171430 у 5B:535723924) -> p2/<g>_extra.paf;
#  2) miniprot белков прочих кандидатов повторно с -G 20k -p 0.3 --outs=0.3 -N 20 (на сервере был --outs по умолчанию 0.99,
#     вторичные выравнивания отбрасывались) -> p2/<g>_win_long.gff.
source "$(dirname "$0")/00_env.sh"; set +e
MP=/opt/micromamba/envs/s16/bin/miniprot
S=$WORK/s16; cd $S/p2
seqkit grep -r -p "TRITD5Bv1G171430_c2" $S/q/queries.fa > $S/q/queries_extra.fa
for w in *_win.fa; do g=${w%_win.fa}
  minimap2 -t 8 -x asm20 -c --secondary=yes -N 50 -p 0.1 $w $S/q/queries_extra.fa > ${g}_extra.paf 2> /dev/null
  $MP -t 8 --gff --trans -G 20k -p 0.3 --outs=0.3 -N 20 $w $S/q/prot_long.fa > ${g}_win_long.gff 2> /dev/null
  echo "$g extra=$(wc -l < ${g}_extra.paf) long=$(grep -c -P '\tmRNA\t' ${g}_win_long.gff)"
done
