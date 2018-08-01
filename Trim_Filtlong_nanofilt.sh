#!/bin/bash

### .fq ###
### -t 22 ###
### -j 22 ###

mkdir -p scaffolds
mkdir -p done_fq
rm -f list
rm -f list_ori
rm -f list_bb
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> list;done
grep 'bb' list > list_bb
grep -v 'bb' list > list_ori
rm list
while read i; do
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i".fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i".fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"q*l1000.fastq
mv "$i"_1.fq.gz ./done_fq/
mv "$i"_2.fq.gz ./done_fq/
done < list_bb
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz && rm "$i"_1.fq.gz
cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
rm "$i"cutLEFT_1.fq
cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz && rm "$i"_2.fq.gz
cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
rm "$i"cutLEFT_2.fq
bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
rm "$i"cutBoth_1.fq
rm "$i"cutBoth_2.fq
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"bbq*l1000.fastq
mv "$i"bb_1.fq.gz ./done_fq/
mv "$i"bb_2.fq.gz ./done_fq/
done < list_ori

