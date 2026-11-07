#!/usr/bin/env bash
# Раздел 17, шаг 4г (локально, WSL; повтор 16_04e для генов части A). Svevo v1: участки s17/q/svevo_windows.tsv (якоря ± 1.1 Мб):
#  5B — из ref11/svevo_5B.fa (скачана в разделе 11); остальные хромосомы — Ensembl REST /sequence/region (целые хромосомы не скачиваются).
#  Имя участка svevo_<chr>_<начало> (16_04f_copies.norm_target переводит координаты в хромосомные).
#  minimap2 -x asm20 -c --secondary=yes -N 50 -p 0.1 запросов s17/q/queries.fa -> s17/svevo/svevo.paf (как 16_04e).
source "$(dirname "$0")/00_env.sh"
O=$WORK/s17/svevo; mkdir -p $O; cd $O
R=https://rest.ensembl.org/sequence/region/triticum_turgidum
tail -n +2 $WORK/s17/q/svevo_windows.tsv | while IFS=$'\t' read grp c s e n l; do
  f=svevo_${c}_${s}.fa
  ok() { [ -s $f ] && [ $(grep -v '>' $f | tr -d '\nACGTNacgtn' | wc -c) -eq 0 ] && [ $(grep -v '>' $f | tr -d '\n' | wc -c) -gt $(( (e-s)/2 )) ]; }
  ok && continue
  if [ $c = 5B ] || [ $c = 5A ]; then samtools faidx $WORK/ref11/svevo_$c.fa svevo_$c:$s-$e | sed "1s/.*/>svevo_${c}_${s}/" > $f
  else for k in 1 2 3 4 5 6; do    # REST иногда отдаёт страницу ошибки HTML — проверка: только ACGTN и длина > половины участка
      curl -s -m 900 "$R/$c:$s..$e:1?content-type=text/x-fasta" | sed "1s/.*/>svevo_${c}_${s}/" > $f; ok && break; sleep $((20*k)); done; fi
  ok || echo "ОШИБКА: участок $f не получен"
  echo "$grp $f $(grep -v '>' $f | tr -d '\n' | wc -c)"
done
cat svevo_*_*.fa > svevo_regions.fa; for f in svevo_*_*.fa; do samtools faidx $f; done
minimap2 -t 8 -x asm20 -c --secondary=yes -N 50 -p 0.1 svevo_regions.fa $WORK/s17/q/queries.fa > svevo.paf 2> mm.log
wc -l svevo.paf
