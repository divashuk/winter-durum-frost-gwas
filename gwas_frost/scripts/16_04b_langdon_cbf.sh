#!/usr/bin/env bash
# Раздел 16, шаг 4б. Копии CBF (DREB1) в Langdon: miniprot белков CBF (CS IPR045277 группы 5, T. monococcum Miller 2006, Svevo IPR045277)
# по окнам Langdon 5A 512–516 и 528–540 Мб, 5B 519–528 Мб. Выход: WSL s16/langdon_cbf/*.gff, *.paf
source "$(dirname "$0")/00_env.sh"
S16=$WORK/s16; O=$S16/langdon_cbf; mkdir -p $O; cd $O
MP=/opt/micromamba/envs/s16/bin/miniprot
P=$WORK/primers11
[ -e $P/langdon_5A.fa.fai ] || samtools faidx $P/langdon_5A.fa
# белки CBF: CS (только IPR045277), Tm, Svevo (если уже скачаны)
cut -f1 $S16/biomart/cs_IPR045277_5ABD.tsv | tail -n +2 > cs_cbf_ids.txt
seqkit grep -n -r -f <(sed 's/^/gene:/;s/$/ /' cs_cbf_ids.txt) $S16/seq/cs_pep.fa | seqkit replace -p '^(\S+) .*gene:(\S+) .*' -r 'CS|$2' > prot.fa
seqkit replace -p '^(\S+) .*binding factor (\d+).*' -r 'Tm|CBF$2|$1' $S16/seq/tm_cbf_prot.fa | seqkit replace -p '^(\S+) TmCBF7.*' -r 'Tm|CBF7|$1' >> prot.fa
[ -s $S16/seq/svevo_pep.fa ] && cut -f1 $S16/biomart/svevo_IPR045277_5A5B.tsv | tail -n +2 > sv_ids.txt && \
  seqkit grep -r -f sv_ids.txt $S16/seq/svevo_pep.fa | seqkit replace -p '^(\S+) .*' -r 'SV|$1' >> prot.fa
grep -c '>' prot.fa
samtools faidx $P/langdon_5A.fa langdon_5A:512000000-516000000 langdon_5A:528000000-540000000 > win.fa
samtools faidx $P/langdon_5B.fa langdon_5B:519000000-528000000 >> win.fa
$MP -t 8 --gff -P CBF -G 3k -p 0.3 --outs=0.3 -N 100 --outn=500 win.fa prot.fa > cbf.gff 2> mp.log
grep -c 'mRNA' cbf.gff
