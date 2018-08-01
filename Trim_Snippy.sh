#!/bin/bash

### Step1: Trim 10bp from both sides and adapters, then discard reads below quality of 20 or length below 2/3 of the original ### 
### Step2: Into Snippy ###
### Change below accordingly ###

### .fq ###
### --cpus 22 ###
### -t 22 ###
### -j 22 ###

mkdir -p done_fq
rm -f list
rm -f list_ori
rm -f list_bb
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> list;done
grep 'bb' list > list_bb
grep -v 'bb' list > list_ori
rm -f list
rm -f list_gbk
for i in *.gbk;do basename -s .gbk "$i" >> list_gbk;done
while read u; do
while read i; do
snippy --cpus 22 --ram 32 --unmapped --outdir ./ref"$u"/"$i"_ref"$u" --ref "$u".gbk --R1 "$i"_1.f*.gz --R2 "$i"_2.f*.gz --bwaopt -t 22
cp ./ref"$u"/"$i"_ref"$u"/snps.vcf ./ref"$u"/"$i"_ref"$u"/"$i"_ref"$u".vcf
cp ./ref"$u"/"$i"_ref"$u"/snps.csv ./ref"$u"/"$i"_ref"$u"/"$i"_ref"$u".csv
mv "$i"_1.f*.gz ./done_fq/
mv "$i"_2.f*.gz ./done_fq/
done < list_bb
done < list_gbk
while read u; do
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz
rm -f "$i"_1.fq.gz
cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
rm -f "$i"cutLEFT_1.fq
cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz
rm -f "$i"_2.fq.gz
cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
rm -f "$i"cutLEFT_2.fq
bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
rm -f "$i"cutBoth_1.fq
rm -f "$i"cutBoth_2.fq
snippy --cpus 22 --ram 32 --unmapped --outdir ./ref"$u"/"$i"bb_ref"$u" --ref "$u".gbk --R1 "$i"bb_1.fq.gz --R2 "$i"bb_2.fq.gz --bwaopt -t 22
cp ./ref"$u"/"$i"bb_ref"$u"/snps.vcf ./ref"$u"/"$i"bb_ref"$u"/"$i"bb_ref"$u".vcf
cp ./ref"$u"/"$i"bb_ref"$u"/snps.csv ./ref"$u"/"$i"bb_ref"$u"/"$i"bb_ref"$u".csv
mv "$i"bb_1.fq.gz ./done_fq/
mv "$i"bb_2.fq.gz ./done_fq/
done < list_ori
done < list_gbk







