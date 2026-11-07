#!/usr/bin/env bash
# Раздел 16, шаг 4к (локально). Второй проход для Svevo v1 (как 16_04i_server_pass2.sh): окна -> miniprot, копии -> fasta.
# Svevo: 5A, 5B — ref11/svevo_5A.fa, svevo_5B.fa (имена svevo_5A/svevo_5B); 7B, 1B — участки ±20 кб (s16/svevo/svevo_<chr>_<начало>.fa).
source "$(dirname "$0")/00_env.sh"
MP=/opt/micromamba/envs/s16/bin/miniprot
S=$WORK/s16; L=$S/pass2/pass2_Svevo.tsv; O=$S/p2; mkdir -p $O; cd $S/svevo
for f in svevo_7B_*.fa svevo_1B_*.fa; do [ -e $f.fai ] || samtools faidx $f; done
loc() { c=$1; s=$2; e=$3; case $c in 5A|5B) echo "$WORK/ref11/svevo_$c.fa svevo_$c:$s-$e $s $e";;
  *) for f in svevo_${c}_*.fa; do o=${f#svevo_${c}_}; o=${o%.fa}; n=$(cut -f2 $f.fai); z=$((o+n-1))
       if [ $e -ge $o ] && [ $s -le $z ]; then a=$(( s>o ? s : o )); b=$(( e<z ? e : z ))   # окно обрезается по участку
       echo "$f ${f%.fa}:$((a-o+1))-$((b-o+1)) $a $b"; return; fi; done;; esac; }
: > $O/Svevo_win.fa; : > $O/Svevo_copies.fa
while IFS=$'\t' read k name c s e st; do
  r=$(loc $c $s $e); [ -z "$r" ] && { echo "нет последовательности для $c:$s-$e"; continue; }
  set -- $r
  if [ $k = win ]; then samtools faidx $1 $2 | sed "1s/.*/>$c:$3-$4/" >> $O/Svevo_win.fa
  else if [ "$st" = "-" ]; then samtools faidx -i $1 $2; else samtools faidx $1 $2; fi | sed "1s/.*/>$name/" >> $O/Svevo_copies.fa; fi
done < $L
$MP -t 8 --gff --trans -G 3k -p 0.3 --outs=0.3 -N 100 --outn=1000 $O/Svevo_win.fa $S/q/prot2.fa > $O/Svevo_win.gff 2> $O/Svevo_mp.log
$MP -t 8 --gff --trans -G 20k -N 5 $O/Svevo_win.fa $S/q/prot_long.fa > $O/Svevo_win_long.gff 2>> $O/Svevo_mp.log
grep -c '>' $O/Svevo_win.fa $O/Svevo_copies.fa; grep -c -P '\tmRNA\t' $O/Svevo_win.gff $O/Svevo_win_long.gff
