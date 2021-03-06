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
spades.py -t 22 --careful --cov-cutoff 10 --pe1-1 "$i"_1.fq.gz --pe1-2 "$i"_2.fq.gz -o "$i"_spades && mkdir -p ./scaffolds/"$i"_spades && cp ./"$i"_spades/scaffolds.fasta ./scaffolds/"$i"_spades/"$i"_spades.fasta && cp ./"$i"_spades/spades.log ./scaffolds/"$i"_spades/"$i"_spades.log && mv "$i"_1.fq.gz ./done_fq/ && mv "$i"_2.fq.gz ./done_fq/
if ([ -f ./"$i"_spades/warnings.log ]);
then
cp ./"$i"_spades/warnings.log ./scaffolds/"$i"_spades/"$i"_spades.warnings && rm -r "$i"_spades
else
rm -r "$i"_spades
fi
done < list_bb
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz && rm "$i"_1.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq && rm "$i"cutLEFT_1.fq && cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz && rm "$i"_2.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq && rm "$i"cutLEFT_2.fq && bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667 && rm "$i"cutBoth_1.fq && rm "$i"cutBoth_2.fq
spades.py -t 22 --careful --cov-cutoff 10 --pe1-1 "$i"bb_1.fq.gz --pe1-2 "$i"bb_2.fq.gz -o "$i"bb_spades && mkdir -p ./scaffolds/"$i"bb_spades && cp ./"$i"bb_spades/scaffolds.fasta ./scaffolds/"$i"bb_spades/"$i"bb_spades.fasta && cp ./"$i"bb_spades/spades.log ./scaffolds/"$i"bb_spades/"$i"bb_spades.log && mv "$i"bb_1.fq.gz ./done_fq/ && mv "$i"bb_2.fq.gz ./done_fq/
if ([ -f ./"$i"bb_spades/warnings.log ]);
then
cp ./"$i"bb_spades/warnings.log ./scaffolds/"$i"bb_spades/"$i"bb_spades.warnings && rm -r "$i"bb_spades
else
rm -r "$i"bb_spades
fi
done < list_ori






