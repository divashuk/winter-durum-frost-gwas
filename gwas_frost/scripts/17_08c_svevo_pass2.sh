#!/usr/bin/env bash
# Раздел 17, шаг 4ж (локально; повтор 16_04j для генов части A). Второй проход для Svevo v1 по участкам s17/svevo/svevo_<chr>_<начало>.fa:
#  окна -> miniprot (параметры 16_04l: -G 20k -p 0.3 --outs=0.3 -N 20), копии -> fasta. Окно обрезается по скачанному участку (как в 16_04j).
source "$(dirname "$0")/00_env.sh"
MP=/opt/micromamba/envs/s16/bin/miniprot
S=$WORK/s17; L=$S/pass2/pass2_Svevo.tsv; O=$S/p2; mkdir -p $O; cd $S/svevo
loc() { c=$1; s=$2; e=$3
  for f in svevo_${c}_*.fa; do o=${f#svevo_${c}_}; o=${o%.fa}; n=$(cut -f2 $f.fai); z=$((o+n-1))
    if [ $e -ge $o ] && [ $s -le $z ]; then a=$(( s>o ? s : o )); b=$(( e<z ? e : z ))
      echo "$f ${f%.fa}:$((a-o+1))-$((b-o+1)) $a $b"; return; fi; done; }
: > $O/Svevo_win.fa; : > $O/Svevo_copies.fa
while IFS=$'\t' read k name c s e st; do
  r=$(loc $c $s $e); [ -z "$r" ] && { echo "нет последовательности для $c:$s-$e"; continue; }
  set -- $r
  if [ $k = win ]; then samtools faidx $1 $2 | sed "1s/.*/>$c:$3-$4/" >> $O/Svevo_win.fa
  else if [ "$st" = "-" ]; then samtools faidx -i $1 $2; else samtools faidx $1 $2; fi | sed "1s/.*/>$name/" >> $O/Svevo_copies.fa; fi
done < $L
: > $O/Svevo_win.gff
$MP -t 8 --gff --trans -G 20k -p 0.3 --outs=0.3 -N 20 $O/Svevo_win.fa $S/q/prot_long.fa > $O/Svevo_win_long.gff 2> $O/Svevo_mp.log
grep -c '>' $O/Svevo_win.fa $O/Svevo_copies.fa; grep -c -P '\tmRNA\t' $O/Svevo_win_long.gff
