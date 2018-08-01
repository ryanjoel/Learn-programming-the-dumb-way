#!/bin/bash

### Step1: Trim 10bp from both sides and adapters, then discard reads below quality of 20 or length below 2/3 of the original ### 
### Step2: Assemble using Unicycler with conservative mode ###
### Change suffix below accordingly ###

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
unicycler -t 22 --mode conservative --min_fasta_length 200 --verbosity 3 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -o "$i"_conservative
mkdir ./scaffolds/"$i"_conservative
cp ./"$i"_conservative/assembly.fasta ./scaffolds/"$i"_conservative/"$i"_conservative.fasta
cp ./"$i"_conservative/assembly.gfa ./scaffolds/"$i"_conservative/"$i"_conservative.gfa
cp ./"$i"_conservative/unicycler.log ./scaffolds/"$i"_conservative/"$i"_conservative.log
rm -r "$i"_conservative
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
unicycler -t 22 --mode conservative --min_fasta_length 200 --verbosity 3 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -o "$i"bb_conservative
mkdir ./scaffolds/"$i"bb_conservative
cp ./"$i"bb_conservative/assembly.fasta ./scaffolds/"$i"bb_conservative/"$i"bb_conservative.fasta
cp ./"$i"bb_conservative/assembly.gfa ./scaffolds/"$i"bb_conservative/"$i"bb_conservative.gfa
cp ./"$i"bb_conservative/unicycler.log ./scaffolds/"$i"bb_conservative/"$i"bb_conservative.log
rm -r "$i"bb_conservative
mv "$i"bb_1.fq.gz ./done_fq/
mv "$i"bb_2.fq.gz ./done_fq/
done < list_ori






