#!/bin/bash

### .fq ###
### -t 22 ###
### -j 22 ###

mkdir -p scaffolds
mkdir -p done_fq
mkdir -p helper_list
rm -f ./helper_list/list
rm -f ./helper_list/list_ori
rm -f ./helper_list/list_bb
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> ./helper_list/list;done
grep 'bb' ./helper_list/list > ./helper_list/list_bb
grep -v 'bb' ./helper_list/list > ./helper_list/list_ori
rm -f ./helper_list/list
while read i; do
if ([ -f "${i:0:6}"*Mde.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i".fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i".fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"q*l1000.fastq && mv "${i:0:6}"*Mde.fastq.gz ./done_fq/ && mv "${i:0:6}"_"$run_info"filton"$i".fastq.gz ./done_fq/
if ([ -f "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq.gz ]);
then
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -l "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq.gz -o "$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative && mkdir -p ./scaffolds/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative && cp ./"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative/assembly.fasta ./scaffolds/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative.fasta && cp ./"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative/assembly.gfa ./scaffolds/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative.gfa && cp ./"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative/unicycler.log ./scaffolds/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative.log && rm -f -r "$i"_"$run_info"filton"$i"q10l1000_hybrid_conservative
mv "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq.gz as it was not found in the current directory"""
fi
else
if ([ -f "${i:0:6}"_*filton*q10l1000.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"_*filton*q10l1000.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -l "${i:0:6}"_"$run_info".fastq.gz -o "$i"_"$run_info"_hybrid_conservative && mkdir -p ./scaffolds/"$i"_"$run_info"_hybrid_conservative && cp ./"$i"_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.fasta && cp ./"$i"_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.gfa && cp ./"$i"_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.log && rm -f -r "$i"_"$run_info"_hybrid_conservative
mv "${i:0:6}"_"$run_info".fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_*filton*q10l1000.fastq.gz as it was not found in the current directory"""
fi
fi
done < ./helper_list/list_bb
###############################################################################################################################################################################
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz && rm "$i"_1.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq && rm "$i"cutLEFT_1.fq && cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz && rm "$i"_2.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq && rm "$i"cutLEFT_2.fq && bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667 && rm "$i"cutBoth_1.fq && rm "$i"cutBoth_2.fq
if ([ -f "${i:0:6}"*Mde.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"bbq*l1000.fastq && mv "${i:0:6}"*Mde.fastq.gz ./done_fq/ && mv "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz ./done_fq/
if ([ -f "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq.gz ]);
then
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq.gz -o "$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative && mkdir -p ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative && cp ./"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative.fasta && cp ./"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative.gfa && cp ./"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative.log && rm -f -r "$i"bb_"$run_info"filton"$i"bbq10l1000_hybrid_conservative
mv "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq.gz as it was not found in the current directory"""
fi
else
if ([ -f "${i:0:6}"_*filton*q10l1000.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"_*filton*q10l1000.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "${i:0:6}"_"$run_info".fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative && mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative && cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta && cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa && cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log && rm -f -r "$i"bb_"$run_info"_hybrid_conservative
mv "${i:0:6}"_"$run_info".fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_*filton*q10l1000.fastq.gz as it was not found in the current directory"""
fi
fi
done < ./helper_list/list_ori
mkdir -p scaffolds
mkdir -p done_fq
mkdir -p helper_list
rm -f ./helper_list/list
rm -f ./helper_list/list_ori
rm -f ./helper_list/list_bb
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> ./helper_list/list;done
grep 'bb' ./helper_list/list > ./helper_list/list_bb
grep -v 'bb' ./helper_list/list > ./helper_list/list_ori
rm -f ./helper_list/list
while read i; do
if ([ -f "${i:0:6}"*Mde.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i".fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i".fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"q*l1000.fastq && mv "${i:0:6}"*Mde.fastq.gz ./done_fq/ && mv "${i:0:6}"_"$run_info"filton"$i".fastq.gz ./done_fq/
if ([ -f "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq.gz ]);
then
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -l "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq.gz -o "$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative && mkdir -p ./scaffolds/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative && cp ./"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative/assembly.fasta ./scaffolds/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative.fasta && cp ./"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative/assembly.gfa ./scaffolds/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative.gfa && cp ./"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative/unicycler.log ./scaffolds/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative.log && rm -f -r "$i"_"$run_info"filton"$i"q9l1000_hybrid_conservative
mv "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq.gz as it was not found in the current directory"""
fi
else
if ([ -f "${i:0:6}"_*filton*q9l1000.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"_*filton*q9l1000.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -l "${i:0:6}"_"$run_info".fastq.gz -o "$i"_"$run_info"_hybrid_conservative && mkdir -p ./scaffolds/"$i"_"$run_info"_hybrid_conservative && cp ./"$i"_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.fasta && cp ./"$i"_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.gfa && cp ./"$i"_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.log && rm -f -r "$i"_"$run_info"_hybrid_conservative
mv "${i:0:6}"_"$run_info".fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_*filton*q9l1000.fastq.gz as it was not found in the current directory"""
fi
fi
done < ./helper_list/list_bb
###############################################################################################################################################################################
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz && rm "$i"_1.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq && rm "$i"cutLEFT_1.fq && cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz && rm "$i"_2.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq && rm "$i"cutLEFT_2.fq && bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667 && rm "$i"cutBoth_1.fq && rm "$i"cutBoth_2.fq
if ([ -f "${i:0:6}"*Mde.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"bbq*l1000.fastq && mv "${i:0:6}"*Mde.fastq.gz ./done_fq/ && mv "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz ./done_fq/
if ([ -f "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq.gz ]);
then
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq.gz -o "$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative && mkdir -p ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative && cp ./"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative.fasta && cp ./"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative.gfa && cp ./"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative.log && rm -f -r "$i"bb_"$run_info"filton"$i"bbq9l1000_hybrid_conservative
mv "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq.gz as it was not found in the current directory"""
fi
else
if ([ -f "${i:0:6}"_*filton*q9l1000.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"_*filton*q9l1000.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "${i:0:6}"_"$run_info".fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative && mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative && cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta && cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa && cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log && rm -f -r "$i"bb_"$run_info"_hybrid_conservative
mv "${i:0:6}"_"$run_info".fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_*filton*q9l1000.fastq.gz as it was not found in the current directory"""
fi
fi
done < ./helper_list/list_ori
mkdir -p scaffolds
mkdir -p done_fq
mkdir -p helper_list
rm -f ./helper_list/list
rm -f ./helper_list/list_ori
rm -f ./helper_list/list_bb
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> ./helper_list/list;done
grep 'bb' ./helper_list/list > ./helper_list/list_bb
grep -v 'bb' ./helper_list/list > ./helper_list/list_ori
rm -f ./helper_list/list
while read i; do
if ([ -f "${i:0:6}"*Mde.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i".fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i".fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"q9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"q10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"q*l1000.fastq && mv "${i:0:6}"*Mde.fastq.gz ./done_fq/ && mv "${i:0:6}"_"$run_info"filton"$i".fastq.gz ./done_fq/
if ([ -f "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq.gz ]);
then
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -l "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq.gz -o "$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative && mkdir -p ./scaffolds/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative && cp ./"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative/assembly.fasta ./scaffolds/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative.fasta && cp ./"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative/assembly.gfa ./scaffolds/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative.gfa && cp ./"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative/unicycler.log ./scaffolds/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative/"$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative.log && rm -f -r "$i"_"$run_info"filton"$i"q8l1000_hybrid_conservative
mv "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_"$run_info"filton"$i"q8l1000.fastq.gz as it was not found in the current directory"""
fi
else
if ([ -f "${i:0:6}"_*filton*q8l1000.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"_*filton*q8l1000.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"_1.fq.gz -2 "$i"_2.fq.gz -l "${i:0:6}"_"$run_info".fastq.gz -o "$i"_"$run_info"_hybrid_conservative && mkdir -p ./scaffolds/"$i"_"$run_info"_hybrid_conservative && cp ./"$i"_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.fasta && cp ./"$i"_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.gfa && cp ./"$i"_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"_"$run_info"_hybrid_conservative/"$i"_"$run_info"_hybrid_conservative.log && rm -f -r "$i"_"$run_info"_hybrid_conservative
mv "${i:0:6}"_"$run_info".fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_*filton*q8l1000.fastq.gz as it was not found in the current directory"""
fi
fi
done < ./helper_list/list_bb
###############################################################################################################################################################################
while read i; do
cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.fq.gz && rm "$i"_1.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq && rm "$i"cutLEFT_1.fq && cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.fq.gz && rm "$i"_2.fq.gz && cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq && rm "$i"cutLEFT_2.fq && bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667 && rm "$i"cutBoth_1.fq && rm "$i"cutBoth_2.fq
if ([ -f "${i:0:6}"*Mde.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"*Mde.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
filtlong -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz --min_length 1000 --trim --split 500 "${i:0:6}"_"$run_info".fastq.gz | pigz > "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz && pigz -d -c "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz | NanoFilt -q 8 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq | NanoFilt -q 9 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq && cat "${i:0:6}"_"$run_info"filton"$i"bbq9l1000.fastq | NanoFilt -q 10 -l 1000 > "${i:0:6}"_"$run_info"filton"$i"bbq10l1000.fastq && pigz "${i:0:6}"_"$run_info"filton"$i"bbq*l1000.fastq && mv "${i:0:6}"*Mde.fastq.gz ./done_fq/ && mv "${i:0:6}"_"$run_info"filton"$i"bb.fastq.gz ./done_fq/
if ([ -f "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq.gz ]);
then
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq.gz -o "$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative && mkdir -p ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative && cp ./"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative.fasta && cp ./"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative.gfa && cp ./"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative/"$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative.log && rm -f -r "$i"bb_"$run_info"filton"$i"bbq8l1000_hybrid_conservative
mv "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_"$run_info"filton"$i"bbq8l1000.fastq.gz as it was not found in the current directory"""
fi
else
if ([ -f "${i:0:6}"_*filton*q8l1000.fastq.gz ]);
then
full_minion_name=$(echo "${i:0:6}"_*filton*q8l1000.fastq.gz) && no_suffix="${full_minion_name%%.*}" && run_info="${no_suffix#*_}"
unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "${i:0:6}"_"$run_info".fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative && mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative && cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta && cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa && cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log && rm -f -r "$i"bb_"$run_info"_hybrid_conservative
mv "${i:0:6}"_"$run_info".fastq.gz ./done_fq/
else
echo """Skip assembling with "${i:0:6}"_*filton*q8l1000.fastq.gz as it was not found in the current directory"""
fi
fi
done < ./helper_list/list_ori
