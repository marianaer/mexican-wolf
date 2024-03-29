#command line eg
# sh MODERN_DSlib_pipeline_single_read_mito_complete_eg_04-06-2015.sh $sample $2 $3 $4 $5

# $1 name of job
#               In order to keep samples reasonably sorted and facilitate coding each sample has
#               a unique prefix (this should be filename for the file filename.fastq.gz) each
#               step in the pipeline then has a suffix
# $2 absolute path to the initial input fastqs
# $3 absolute path to the output directory
#               IMPORTANT:  Do not include a "/" after the last directory it will break this script
# $4 path and filename of the reference genome ie. /mnt/tank/species/Dog/CanFam2.fa
#               IMPORTANT:  The reference genome must be correctly indexed with bwa index and samtools faidx
# $5 bwa mismatch value:  allows easy comparisons between different mismatch values (0.04)
# $6 directory of the sample
# $7 name of the sample

# ### paths the executable files used in the script
Picard=/raid6-18tb/software/picard-1.111/AddOrReplaceReadGroups.jar
fastq_screen=/raid6/utaron/software/fastq_screen_v0.5.0/fastq_screen
fastqc=/raid6/utaron/software/FastQC/fastqc

### Variables - set these to something appropriate for the sample
MAPQ=30                         # Minimum read mapping quality for bwa
Mismatch=$5                     # Mismatch threshold
                                # 0.04 (default value for bwa) is 5 mismatches in 100 bp
                                # 0.01 is ~6 mismatches in 100 bp
                                # bwa aln prints the number of mismatches allowed at the beginning of the run
max_insert=1000                 # for PE mapping: bwa sampe uses this information for mapping unmerged PE reads
threads=5                     # no. of computer cores used by bwa and samtools. 20 = OK, >20 = ask people first!; max 2 on your own machine


###### MAIN SCRIPT ######

### Create ouput diretory
mkdir $3
mkdir $3/$6

# ### trim adapter seqs using cutadapt
# # open programme on computing server: cutadapt
# # -b = adapter sequence
# # -m = min read length; shorter sequences will not be incoporated
# # -o = output file .fastq
cutadapt-1.12 -b AGATCGGAAGAGCACACGTC -m 30 -o $3/$6/$7.cutadapt.fastq $2/$6/"$7"_R1.fastq


#
#
# ### zip initial fastq files
gzip $2/$6/"$7"_R1.fastq

#echo $2/$6/"$7"_R1.fastq

# ### zip trim.fastq files
# # compress cutadapt.fastq as it is a huge file
gzip $3/$6/"$7".cutadapt.fastq


#
# ### check reads in FastQScreen
# # Fastq_Screen version 5 database = nuclear genomes
fastq_screen-v0.4.4 -threads $threads -subset 1000000 -outdir $3/$6 -conf /raid6/mescobar/wolves/config/fastq_screen_nuclear.conf  -aligner bowtie2 $3/$6/"$7".cutadapt.fastq.gz

### Read alignment and mapping for singe end reads: bwa aln samse, samtools view and sort
bwa-v0.7.8 aln -n $Mismatch -t $threads $4 $3/$6/"$7".cutadapt.fastq.gz | bwa-v0.7.8 samse $4 - $3/$6/"$7".cutadapt.fastq.gz | samtools-v0.1.19 view -Su -q $MAPQ -@ $threads - | samtools-v0.1.19 sort -@ $threads -m 8G - $3/$6/"$7".nucl.cutadapt_sort

## Picard MarkDuplicates. Remove dupliates from the .bam file
java -jar /raid6-18tb/software/picard-1.111/MarkDuplicates.jar I=$3/$6/"$7".nucl.cutadapt_sort.bam O=$3/$6/"$7".nucl.cutadapt_rmdup_unsort.bam M=marked_dup_metrics.txt VALIDATION_STRINGENCY=SILENT

### samtools sort: sort the resulting .bam file
samtools-v0.1.19 sort -@ $threads $3/$6/"$7".nucl.cutadapt_rmdup_unsort.bam $3/$6/"$7".nucl.cutadapt_rmdup_sort

### samtools index: index the resulting .bam file
samtools-v0.1.19 index $3/$6/"$7".nucl.cutadapt_rmdup_sort.bam

### picardtools adding read groups
java -jar -XX:ParallelGCThreads=$threads $Picard INPUT=$3/$6/"$7".nucl.cutadapt_rmdup_sort.bam OUTPUT=$3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam RGID="$7" RGLB="$7"_rglb RGPL=Illumina RGPU="$7"_rgpu RGSM="$7"_rgsm VALIDATION_STRINGENCY=SILENT

### index picardtools output
samtools-v0.1.19 index $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam

#### report number of mapped and unmapped nuclear reads
samtools-v0.1.19 idxstats $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam | awk 'BEGIN {a=0;b=0} {a += $3; b+=$4 } END{print " nuclear mapped $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam " a " nulcear unmapped $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam " b}'

### report average nuclear read depth and total mapped bp
#samtools-v0.1.19 depth $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam | awk '{sum+=$3;cnt++}END{print " nulcear read depth $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam  " sum/cnt " total nulcear mapped bp $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam " sum}'

### run mapDamge on the mapped reads
mapDamage --no-stats -i $3/$6/"$7".nucl.cutadapt_rmdup_sort_rg.bam -r $4 -d $3/$6/"$7"_nucl

echo $1 done
