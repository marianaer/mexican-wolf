#!/bin/bash

cd $2
samples=(Jal545 Jal547 Jal548 Jal549 Jal550)

for direct in *; do
mkdir $3/$direct
done

# Loop through each sample in each directory
for sample in "${samples[@]}"; do
	echo  $sample
	for direct in *; do
           
        	if [[ $direct =~ .*"$sample".* ]] && [[ $direct =~ .*_14.* || $direct =~ .*_30.* ]]; then  # If the reads are paired end
             
             		zcat $direct/$sample*R1*.fastq.gz > /$2/$direct/"$sample"_R1.fastq
             		zcat $direct/$sample*R2*.fastq.gz > /$2/$direct/"$sample"_R2.fastq
	     		/raid6/mescobar/wolves/dstest.sh $1 $2 $3 $4 $5 $direct $sample & # map with paired end script
	


		elif [[ $direct =~ .*"$sample".* ]]; then # if they are single end instead

           		zcat $direct/$sample*R1*.fastq.gz > /$2/$direct/"$sample"_R1.fastq
           		/raid6/mescobar/wolves/sstest.sh $1 $2 $3 $4 $5 $direct $sample & # map with single end sript
 			
	
        	fi
	done
done


# Merge the resulting mapped bam files...

cd $3

for sample in "${samples[@]}"; do  
       mkdir $3/$sample

        arr=(`find /raid6/mescobar/wolves/outputs -type f -name "$sample*rg.bam" | sort -r | head`)  # find the .bam files for each of the samples
	#echo "${array[*]}"
	
        # Merging each sample's bam files
        java -jar /raid6-18tb/software/picard-1.111/MergeSamFiles.jar   I="${arr[0]}"  I="${arr[1]}" I="${arr[2]}"  O=$3/$sample/"$sample"_se_pe_mergedfghfghgfhfghgfh.bam VALIDATION_STRINGENCY=SILENT


        # Sorting the resuting merged bam file
        samtools-v0.1.19 sort -@10 $3/$sample/"$sample"_se_pe_merged.bam $3/$sample/"$sample"_se_pe_merged_sort

        # Indexing the resulting merged bam file
        samtools-v0.1.19 index $3/$sample/"$sample"_se_pe_merged_sort.bam

done



