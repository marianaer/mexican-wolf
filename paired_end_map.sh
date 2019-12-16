#########################################
###### PAIRED END LIBRARY MAPPINGG ######
#########################################

###### Overview ######
# Pipeline for read merging and bwa mapping to nulcear reference
# Illumina reads
# paired-end
# dsDNA libraries
# historic samples
# Only Merged PE reads are mapped (if needed unmerged reads can be included)

##### Command Line arguments #####

## command line eg
# sh MODERN_DSlib_pipeline_eg_04-06-2015.sh $1 $2 $3 $4 $5

# $1 core name of input file
#               In order to keep samples reasonably sorted and facilitate coding each sample has
#               a unique prefix (this should be filename for the file filename.fastq.gz) each
#               step in the pipeline then has a suffix
# $2 path to the initial input fastqs
# $3 path to the output directory
#               IMPORTANT:  Do not include a "/" after the last directory it will break this script
# $4 path and filename of the reference genome ie. /mnt/tank/species/Dog/CanFam2.fa
#               IMPORTANT:  The reference genome must be correctly indexed with bwa index and samtools faidx
# $5 bwa mismatch value:  allows easy comparisons between different mismatch values

### paths the executable files used in the script
Picard=/raid6-18tb/software/picard-1.111/AddOrReplaceReadGroups.jar
fastq_screen=/raid6/utaron/software/fastq_screen_v0.5.0/fastq_screen
fastqc=/raid6/utaron/software/FastQC/fastqc
flash=/raid6/mescobar/software/FLASH-1.2.11/flash

### Variables - set these to something appropriate for the sample
MAPQ=30               		# Minimum read mapping quality for bwa
Mismatch=$5  			# Mismatch threshold
                                # 0.04 (default value for bwa) is 5 mismatches in 100 bp
                                # 0.01 is ~6 mismatches in 100 bp
                                # bwa aln prints the number of mismatches allowed at the beginning of the run
merge_overlap=15		# minimum overlap for SeqPrep merging
merge_MinLen=25			# minimum accepted read length (SeqPrep)
merge_qual=13			# minimum merge quality (SeqPrep)
max_insert=1000			# for PE mapping: bwa sampe uses this information for mapping unmerged PE reads
threads=10                      # no. of computer cores used by bwa and samtools. The local machine has 4 cores, so use less than 4 if you want to work in parallel while the script is running.


###### MAIN SCRIPT ######

### Create output directory
mkdir $3

### concatenate all R1 and R2 fastqs for the sample (assumes files are gzipped)
zcat $2/$1*R1*.fastq.gz > $2/$1_R1.fastq
zcat $2/$1*R2*.fastq.gz > $2/$1_R2.fastq

# ### inactivation
# ##trim adapter seqs and merge overlapping R1 and R2

# ### Adapter trimming
cutadapt-1.12 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT -o $2/$1_out_1.cutadapt.fastq -p $2/$1_out_2.cutadapt.fastq  $2/$1_R1.fastq $2/$1_R2.fastq -j 5

# # Read merging
$flash $2/$1_out_1.cutadapt.fastq $2/$1_out_2.cutadapt.fastq -m 15 -o $2/$1_out_flash -M 75 -d $3


# ### inactivation
# ##zip initial fastq files and resulting merged file
gzip $2/$1_R1.fastq
gzip $2/$1_R2.fastq
gzip $2/$1_out_flash*

########################################
### Is is very likely that most of the reads are overlapping and will be merged. That is why the downstream analysis focusses only on the merged reads.
########################################

# ### inspect reads after adapter removal with FastQC
# # --noextract = outputfile will not be uncompressed
# # -o = output directory
#fastqc $2/$1_out_flash_extendedFrags.fastq.gz --noextract -o $3
#
# ### check reads in FastQScreen
# # Fastq_Screen version 5 database = nuclear genomes
fastq_screen-v0.4.4 -threads $threads -subset 1000000 -outdir $3 -conf /raid6/mescobar/wolves/config/fastq_screen_nuclear.conf -aligner bowtie2 $2/$1_out_flash.extendedFrags.fastq.gz

### Read alignment and mapping for merged reads: bwa aln samse, samtools view and sort
## using bwa v 0.6.2 since v 0.7.5a gave some errors
#bwa-v0.7.8 aln -n $Mismatch -t $threads $4 $2/$1_out_flash.extendedFrags.fastq.gz | bwa-v0.7.8 samse
$4 - $2/$1_out_flash.extendedFrags.fastq.gz | samtools-v0.1.19 view -Su -q $MAPQ -@ $threads - | samtools-v0.1.19 sort -@ $threads -m 8G - $3/$1_nucl_pe_merged_reads_sort

# Picard MarkDuplicates
java -jar picard.jar MarkDuplicates I=$3/$1_nucl_pe_merged_reads_sort.bam O=$3/$1_nucl_pe_merged_reads_rmdup_unsort.bam M=marked_dup_metrics.txt

### samtools sort: sort the resulting .bam file
samtools-v0.1.19 sort -@ $threads $3/$1_nucl_pe_merged_reads_rmdup_unsort.bam $3/$1_nucl_pe_merged_reads_rmdup_sort

## samtools index: index the resulting .bam file
samtools-v0.1.19 index $3/$1_nucl_pe_merged_reads_rmdup_sort.bam

## picardtools adding read groups
java -jar -Xmx10240M $Picard INPUT=$3/$1_nucl_pe_merged_reads_rmdup_sort.bam OUTPUT=$3/$1_nucl_pe_merged_reads_rmdup_sort_rg.bam RGID=$1 RGLB=$1_rglb RGPL=Illumina RGPU=$1_rgpu RGSM=$1_rgsm VALIDATION_STRINGENCY=SILENT

## index picardtools output
samtools-v0.1.19 index $3/$1_nucl_pe_merged_reads_rmdup_sort_rg.bam

# report number of mapped and unmapped nuclear reads
samtools-v0.1.19 idxstats $3/$1_nucl_pe_merged_reads_rmdup_sort_rg.bam | awk 'BEGIN {a=0;b=0} {a += $3;b+=$4 } END{print " nucl mapped " a " nucl unmapped " b}'

# report average nuclear read depth and total mapped bp
samtools-v0.1.19 depth $3/$1_nucl_pe_merged_reads_rmdup_sort_rg.bam | awk '{sum+=$3;cnt++}END{print " nucl read depth " sum/cnt " total nucl mapped bp " sum}'

### run mapDamge on the mapped reads
mapDamage --no-stats -i $3/$1_nucl_pe_merged_reads_rmdup_sort_rg.bam -r $4 -d $3/$1_nucl

### Bam merging: peDNA and seDNA
#Before variant calling... Picard reference preparation should be done. Also, change the reference extension (.fna) to .fa , otherwise, Picard won't work

java -jar /raid6-18tb/software/picard-1.111/CreateSequenceDictionary.jar R= $4 O= $4.dict

### BAM MERGING OF BOTH PAIRED END AND SINGLE END READS

java -jar /raid6-18tb/software/picard-1.111/MergeSamFiles.jar   I=$3/$1_nucl_pe_merged_reads_rmdup_sort_rg.bam  I=.$3/$1_S16.nucl.cutadapt_se_rmdup_sort_rg.bam  O=$3/$1_nucl_pe_se_merged.bam VALIDATION_STRINGENCY=SILENT


## Samtools sort and indexing of resulting merged file:

### samtools sort: sort the resulting .bam file
samtools-v0.1.19 sort -@ 10 $3/$1_nucl_pe_se_merged.bam $3/$1_nucl_pe_se_merged_sort

## samtools index: index the resulting .bam filfe
samtools-v0.1.19 index $3/$1_nucl_pe_se_merged_sort.bam

## Variant calling:
java -jar -Xmx10240M /raid6-18tb/software/GATK-2.8-1/GenomeAnalysisTK.jar -I $3/$1_nucl_pe_se_merged_sort.bam -o $3/$1_nucl_pe_se_merged_sort.bam_variants.vcf -R $4 -T UnifiedGenotyper --genotype_likelihoods_model BOTH --sample_ploidy 1 --output_mode EMIT_ALL_SITES
