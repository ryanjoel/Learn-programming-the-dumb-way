#!/bin/bash

# Trim 10bp from both sides and adapters, then discard reads below quality of 20 or length below 2/3 of the original
# Assemble using Unicycler with conservative mode
# Change suffix below accordingly

# .fq
# -j 22
# -t 22

# Preparing list of samples

rm -f -f list
rm -f list_final
for i in *_1.fq.gz;do basename -s _1.fq.gz "$i" >> list;done
sed 's/bb//g' list | sort | uniq > list_final
rm -f list
mv list_final list
mkdir -p scaffolds
mkdir -p done_fq

# Trim and Assembly

while read i;do
    if ([ -f "$i"_*Mdeq10l1000.fastq.gz ]);
    then
        echo " Starting to assemble strain "$i" using hybrid strategy provided by Unicycler with MinION reads (Quality > Q10 and Length > 1000bp) "   # Hybrid assembly with q10l1000 MinION reads
        full_minion_name=$(echo "$i"_*Mdeq10l1000.fastq.gz)
        no_suffix="${full_minion_name%%.*}"
        run_info="${no_suffix#*_}"  
        if ([ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log ]);
        then
            echo " "$i"bb_"$run_info"_hybrid_conservative already exists "
            mv "$i"_*Mdeq10l1000.fastq.gz ./done_fq/
        else
            if ([ -f "$i"bb_1.f*.gz ] && [ -f "$i"bb_2.f*.gz ]);   # If trimmed reads exist
            then
                echo " Trimmed Illumina reads of "$i" already exist, proceeding to assembly with conservative mode "
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mdeq10l1000.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative
                mv "$i"_*Mdeq10l1000.fastq.gz ./done_fq/
            elif ([ -f "$i"_1.f*.gz ] && [ -f "$i"_2.f*.gz ]);   # If trimmed reads do not exit but their original reads exist 
            then
                echo " Trimmed Illumina reads of "$i" do not exist, starting with the original now "       
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
                rm -f "$i"cutLEFT_1.fq
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
                rm -f "$i"cutLEFT_2.fq
                bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
                rm -f "$i"cutBoth_1.fq
                rm -f "$i"cutBoth_2.fq
                rm -f "$i"_1.f*.gz
                rm -f "$i"_2.f*.gz
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mdeq10l1000.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative
                mv "$i"_*Mdeq10l1000.fastq.gz ./done_fq/
            else
                echo " Please provide the paired Illumina reads of "$i" "
            fi
        fi
    else
        echo " MinION Mdeq10l1000 reads of strain "$i" can not be found "  
    fi
done < list
while read i;do
    if ([ -f "$i"_*Mdeq9l1000.fastq.gz ]);
    then
        echo " Starting to assemble strain "$i" using hybrid strategy provided by Unicycler with MinION reads (Quality > Q9 and Length > 1000bp) "   # Hybrid assembly with q9l1000 MinION reads
        full_minion_name=$(echo "$i"_*Mdeq9l1000.fastq.gz)
        no_suffix="${full_minion_name%%.*}"
        run_info="${no_suffix#*_}"  
        if ([ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log ]);
        then
            echo " "$i"bb_"$run_info"_hybrid_conservative already exists "
            mv "$i"_*Mdeq9l1000.fastq.gz ./done_fq/
        else
            if ([ -f "$i"bb_1.f*.gz ] && [ -f "$i"bb_2.f*.gz ]);   # If trimmed reads exist
            then
                echo " Trimmed Illumina reads of "$i" already exist, proceeding to assembly with conservative mode "
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mdeq9l1000.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative
                mv "$i"_*Mdeq9l1000.fastq.gz ./done_fq/
            elif ([ -f "$i"_1.f*.gz ] && [ -f "$i"_2.f*.gz ]);   # If trimmed reads do not exit but their original reads exist 
            then
                echo " Trimmed Illumina reads of "$i" do not exist, starting with the original now "       
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
                rm -f "$i"cutLEFT_1.fq
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
                rm -f "$i"cutLEFT_2.fq
                bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
                rm -f "$i"cutBoth_1.fq
                rm -f "$i"cutBoth_2.fq
                rm -f "$i"_1.f*.gz
                rm -f "$i"_2.f*.gz
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mdeq9l1000.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative
                mv "$i"_*Mdeq9l1000.fastq.gz ./done_fq/
            else
                echo " Please provide the paired Illumina reads of "$i" "
            fi
        fi
    else
        echo " MinION Mdeq9l1000 reads of strain "$i" can not be found "  
    fi
done < list
while read i;do
    if ([ -f "$i"_*Mdeq8l1000.fastq.gz ]);
    then
        echo " Starting to assemble strain "$i" using hybrid strategy provided by Unicycler with MinION reads (Quality > Q8 and Length > 1000bp) "   # Hybrid assembly with q8l1000 MinION reads
        full_minion_name=$(echo "$i"_*Mdeq8l1000.fastq.gz)
        no_suffix="${full_minion_name%%.*}"
        run_info="${no_suffix#*_}"  
        if ([ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log ]);
        then
            echo " "$i"bb_"$run_info"_hybrid_conservative already exists "
            mv "$i"_*Mdeq8l1000.fastq.gz ./done_fq/
        else
            if ([ -f "$i"bb_1.f*.gz ] && [ -f "$i"bb_2.f*.gz ]);   # If trimmed reads exist
            then
                echo " Trimmed Illumina reads of "$i" already exist, proceeding to assembly with conservative mode "
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mdeq8l1000.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative 
                mv "$i"_*Mdeq8l1000.fastq.gz ./done_fq/
            elif ([ -f "$i"_1.f*.gz ] && [ -f "$i"_2.f*.gz ]);   # If trimmed reads do not exit but their original reads exist 
            then
                echo " Trimmed Illumina reads of "$i" do not exist, starting with the original now "       
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
                rm -f "$i"cutLEFT_1.fq
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
                rm -f "$i"cutLEFT_2.fq
                bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
                rm -f "$i"cutBoth_1.fq
                rm -f "$i"cutBoth_2.fq
                rm -f "$i"_1.f*.gz
                rm -f "$i"_2.f*.gz
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mdeq8l1000.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative                            
                mv "$i"_*Mdeq8l1000.fastq.gz ./done_fq/
            else
                echo " Please provide the paired Illumina reads of "$i" "
            fi
        fi
    else
        echo " MinION Mdeq8l1000 reads of strain "$i" can not be found "  
    fi
done < list
while read i;do
    if ([ -f "$i"_*Mde.fastq.gz ]);
    then
        echo " Starting to assemble strain "$i" using hybrid strategy provided by Unicycler with demultiplexed MinION reads "   # Hybrid assembly with original  Mde MinION reads
        full_minion_name=$(echo "$i"_*Mde.fastq.gz)
        no_suffix="${full_minion_name%%.*}"
        run_info="${no_suffix#*_}"  
        if ([ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa ] && [ -f ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log ]);
        then
            echo " "$i"bb_"$run_info"_hybrid_conservative already exists "
            mv "$i"_*Mde.fastq.gz ./done_fq/
        else
            if ([ -f "$i"bb_1.f*.gz ] && [ -f "$i"bb_2.f*.gz ]);   # If trimmed reads exist
            then
                echo " Trimmed Illumina reads of "$i" already exist, proceeding to assembly with conservative mode "
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mde.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative
                mv "$i"_*Mde.fastq.gz ./done_fq/
            elif ([ -f "$i"_1.f*.gz ] && [ -f "$i"_2.f*.gz ]);   # If trimmed reads do not exit but their original reads exist 
            then
                echo " Trimmed Illumina reads of "$i" do not exist, starting with the original now "       
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_1.fq "$i"_1.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_1.fq "$i"cutLEFT_1.fq
                rm -f "$i"cutLEFT_1.fq
                cutadapt -j 22 -u 10 -o "$i"cutLEFT_2.fq "$i"_2.f*.gz
                cutadapt -j 22 -u -10 -o "$i"cutBoth_2.fq "$i"cutLEFT_2.fq
                rm -f "$i"cutLEFT_2.fq
                bbduk.sh in1="$i"cutBoth_1.fq in2="$i"cutBoth_2.fq out1="$i"bb_1.fq.gz out2="$i"bb_2.fq.gz ref=/home/zong/programs_and_scripts/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 maq=20 mlf=0.6666666667
                rm -f "$i"cutBoth_1.fq
                rm -f "$i"cutBoth_2.fq
                rm -f "$i"_1.f*.gz
                rm -f "$i"_2.f*.gz
                unicycler -t 22 --mode conservative --min_fasta_length 200 -1 "$i"bb_1.fq.gz -2 "$i"bb_2.fq.gz -l "$i"_*Mde.fastq.gz -o "$i"bb_"$run_info"_hybrid_conservative
                mkdir -p ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.fasta ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.fasta
                cp ./"$i"bb_"$run_info"_hybrid_conservative/assembly.gfa ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.gfa
                cp ./"$i"bb_"$run_info"_hybrid_conservative/unicycler.log ./scaffolds/"$i"bb_"$run_info"_hybrid_conservative/"$i"bb_"$run_info"_hybrid_conservative.log
                rm -f -r "$i"bb_"$run_info"_hybrid_conservative                
                mv "$i"_*Mde.fastq.gz ./done_fq/
            else
                echo " Please provide the paired Illumina reads of "$i" "
            fi
        fi
    else
        echo " MinION Mde reads of strain "$i" can not be found "  
    fi
done < list
echo " Now the program has reached its end "
echo " Thanks "
