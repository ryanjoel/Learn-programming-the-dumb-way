#!/bin/bash

mkdir -p done_fq
rm -f list
rm -f list_ori
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> list;done
grep -v 'bb' list > list_ori
rm list
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz && mv "$i"_1.fq.gz ./done_fq/
cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq && rm "$i"cutLEFT_1.fq
cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz && mv "$i"_2.fq.gz ./done_fq/
cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq && rm "$i"cutLEFT_2.fq
bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667 && rm "$i"cutBoth_1.fq && rm "$i"cutBoth_2.fq
counter=1
while [ $counter -le 10 ]; do
echo """
#####################################
     Trimming iteration "$counter" starts     
#####################################
"""
bbduk.sh in1="$i"bb_1.fq.gz in2="$i"bb_2.fq.gz out1="$i"bbbb_1.fq.gz out2="$i"bbbb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667 |& tee "$i"_trim.log && rm "$i"bb_1.fq.gz && rm "$i"bb_2.fq.gz && mv "$i"bbbb_1.fq.gz "$i"bb_1.fq.gz && mv "$i"bbbb_2.fq.gz "$i"bb_2.fq.gz
tail -n 11 "$i"_trim.log | head -n 1 > "$i"_last_round.txt && sed 's/ /_/g' "$i"_last_round.txt | tsv2csv.py | csvcut -c 4 > "$i"_last_round_nt_raw.txt && u=$(cat "$i"_last_round_nt_raw.txt) && echo "${u%%_*}" > "$i"_last_round_nt && rm "$i"_last_round.txt && rm "$i"_last_round_nt_raw.txt
tail -n 5 "$i"_trim.log | head -n 1 > "$i"_this_round.txt && sed 's/ /_/g' "$i"_this_round.txt | tsv2csv.py | csvcut -c 3 > "$i"_this_round_nt_raw.txt && u=$(cat "$i"_this_round_nt_raw.txt) && echo "${u%%_*}" > "$i"_this_round_nt && rm "$i"_this_round.txt && rm "$i"_this_round_nt_raw.txt
if ([ $(cat "$i"_this_round_nt) == $(cat "$i"_last_round_nt) ]); then
echo """
####################################################################
     Trimming finished as no change was found after iteration "$counter"     
####################################################################
"""
rm "$i"_last_round_nt
rm "$i"_this_round_nt
rm "$i"_trim.log
break
else
((counter++))
rm "$i"_last_round_nt
rm "$i"_this_round_nt
rm "$i"_trim.log
fi
done
done < list_ori
