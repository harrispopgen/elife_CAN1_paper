#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -N positionFreqs
#$ -l mfree=4G


## takes suffix
### bash positionFreqs.sh  sample_prefix   fastq_PATH  output_PATH  REF  sample_suffix

## used in centos 6
#module load python/2.7.3 samtools/1.3.1 bowtie2/2.2.3  pear/0.9.5  cutadapt/1.8.3  trim_galore/0.4.1  pysam/0.8.4 pysamstats/0.24.2

## update with centos 7
## modules are installed inside anaconda

module load  samtools/1.3.1 bowtie2/2.2.3  pear/0.9.11  fastx-toolkit/0.0.14


SAMPLE=$1
FQ_PATH=$2
OUT_PATH=$3
REF=$4
SUFFIX=lang
Sample=${SAMPLE}_${SUFFIX}
Sample_suffix=$5


## unzip (first time running this script)
gunzip ${FQ_PATH}/${SAMPLE}_R1_${Sample_suffix}.fastq.gz
gunzip ${FQ_PATH}/${SAMPLE}_R2_${Sample_suffix}.fastq.gz



trim_galore --nextera ${FQ_PATH}/${SAMPLE}_R2_${Sample_suffix}.fastq \
     ${FQ_PATH}/${SAMPLE}_R1_${Sample_suffix}.fastq --length 20 --paired -o  $OUT_PATH



pear -j 16 -r ${OUT_PATH}/${SAMPLE}_R2_${Sample_suffix}_val_1.fq -f ${OUT_PATH}/${SAMPLE}_R1_${Sample_suffix}_val_2.fq \
     -o ${OUT_PATH}/${SAMPLE}.trimmed.pear > ${OUT_PATH}/${SAMPLE}.pear.log

cat ${OUT_PATH}/${SAMPLE}.trimmed.pear.unassembled.forward.fastq  \
     ${OUT_PATH}/${SAMPLE}.trimmed.pear.unassembled.reverse.fastq \
     ${OUT_PATH}/${SAMPLE}.trimmed.pear.assembled.fastq > ${OUT_PATH}/${SAMPLE}.trimmed.pear.all.fastq


## will filter fastq reads quality of raw reads
## https://wikis.utexas.edu/display/bioiteam/FASTQ+Manipulation+Tools
fastq_quality_filter -q 20 -p 94 -i ${OUT_PATH}/${SAMPLE}.trimmed.pear.all.fastq    -Q 33 -o  ${OUT_PATH}/${SAMPLE}.trimmed.pear.all.filter.fastq


bowtie2 -p 16 -t --end-to-end --no-unal  -x ../data/CAN1_Lang \
     -U ${OUT_PATH}/${SAMPLE}.trimmed.pear.all.filter.fastq -S ${OUT_PATH}/${Sample}.trimmed.pear.all.sam 2> ${OUT_PATH}/${Sample}.bowtie.log



samtools view -bS ${OUT_PATH}/${Sample}.trimmed.pear.all.sam \
     | samtools sort > ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam

samtools index -b ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam \
     > ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam.bai

samtools depth -d 10000000000 -a ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam \
     > ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.readdepth.txt



## only keep uniquely mapped  reads (for bowtie bam file, MAPQ=0 if it maps to multiple places.  bwa will be more complicated,
## see: /net/harris/vol1/project/yeast_1002/call_var_raw/TWN_CHN/keep_unique_reads.sh)
## and use MAPQ cutoff of 20 (only include uniquely mapped reads + map quality cutoff)

# use MQ 20 and 40 separately (for indels and SNPs)

samtools view -b -q 20  ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam >  ${OUT_PATH}/${Sample}.filter20.bam

samtools index ${OUT_PATH}/${Sample}.filter20.bam


samtools view -b -q 40  ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam >  ${OUT_PATH}/${Sample}.filter40.bam

samtools index ${OUT_PATH}/${Sample}.filter40.bam



## output original stats for comparison
pysamstats -t variation -d -D 10000000 -f $REF \
    ${OUT_PATH}/${Sample}.trimmed.pear.all.sort.bam  > ${OUT_PATH}/${Sample}.all.stats.txt


pysamstats -t variation -d -D 10000000 -f $REF \
    ${OUT_PATH}/${Sample}.filter20.bam  > ${OUT_PATH}/${Sample}.filter20.all.stats.txt

pysamstats -t variation -d -D 10000000 -f $REF \
    ${OUT_PATH}/${Sample}.filter40.bam  > ${OUT_PATH}/${Sample}.filter40.all.stats.txt




#remove no-longer-needed files.

rm ${OUT_PATH}/${SAMPLE}.trimmed.pear.unassembled.forward.fastq
rm ${OUT_PATH}/${SAMPLE}.trimmed.pear.unassembled.reverse.fastq

rm ${OUT_PATH}/${SAMPLE}_R1_${Sample_suffix}_val_2.fq
rm ${OUT_PATH}/${SAMPLE}_R2_${Sample_suffix}_val_1.fq
