#!/bin/bash

mkdir -p done_fq
rm -f list
rm -f list_ori
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> list;done
grep -v 'bb' list > list_ori
rm list
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz
mv "$i"_1.fq.gz ./done_fq/
cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
rm "$i"cutLEFT_1.fq
cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz
mv "$i"_2.fq.gz ./done_fq/
cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
rm "$i"cutLEFT_2.fq
bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
rm "$i"cutBoth_1.fq
rm "$i"cutBoth_2.fq
done < list_ori
