#!/bin/bash

### .fq ###
### -t 22 ###
### -j 22 ###

mkdir -p log
mkdir -p done
rm -f list
rm -f list_ori
rm -f list_bb
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> list;done
grep 'bb' list > list_bb
grep -v 'bb' list > list_ori
rm list
while read i; do
bbmap.sh ref="$i"_conservative.fasta nodisk in1="$i"_1.fq.gz in2="$i"_2.fq.gz unpigz pigz basecov=basecov.txt overwrite |& tee "$i".log
rm basecov.txt
tail -n 6 "$i".log | head -n 1 | tsv2csv.py | csvcut -c 2 > cov.csv
mv "$i".log ./log/
echo "$i"_conservative.fasta,"$(cat cov.csv)" > "$i"_cov.csv
rm cov.csv
mv "$i"_conservative.fasta ./done/
mv "$i"_1.fq.gz ./done/
mv "$i"_2.fq.gz ./done/
done < list_bb
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz
rm "$i"_1.fq.gz
cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
rm "$i"cutLEFT_1.fq
cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz
rm "$i"_2.fq.gz
cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
rm "$i"cutLEFT_2.fq
bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
rm "$i"cutBoth_1.fq
rm "$i"cutBoth_2.fq
bbmap.sh ref="$i"bb_conservative.fasta nodisk in1="$i"bb_1.fq.gz in2="$i"bb_2.fq.gz unpigz pigz basecov=basecov.txt overwrite |& tee "$i"bb.log
rm basecov.txt
tail -n 6 "$i"bb.log | head -n 1 | tsv2csv.py | csvcut -c 2 > cov.csv
mv "$i"bb.log ./log/
echo "$i"bb_conservative.fasta,"$(cat cov.csv)" > "$i"bb_cov.csv
rm cov.csv
mv "$i"bb_conservative.fasta ./done/
mv "$i"bb_1.fq.gz ./done/
mv "$i"bb_2.fq.gz ./done/
done < list_ori
echo ',Average_coverage' > header.csv
cat header.csv *_cov.csv > Average_coverage.csv
rm header.csv
rm *_cov.csv


