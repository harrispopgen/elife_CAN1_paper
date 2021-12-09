#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -N positionFreqs
#$ -l mfree=5G
#$ -l h_rt=5:0:0

### this script use sample fasta as the reference seq 

### modifed on 1.30.20, will add "r2" as suffix to distinguish

## takes suffix
### bash positionFreqs.sh  sample_prefix   fastq_PATH  output_PATH  REF_prefix  sample_suffix

## old centos 6
#module load python/2.7.3 samtools/1.3.1 bowtie2/2.2.3  pear/0.9.5  cutadapt/1.8.3  trim_galore/0.4.1  pysam/0.8.4 pysamstats/0.24.2


module load  samtools/1.3.1 bowtie2/2.2.3  pear/0.9.11

SAMPLE=$1
FQ_PATH=$2
OUT_PATH=$3
REF_dir="../data/fasta"
REF_bowtie="${REF_dir}/${4}_CAN1_Lang"

REF=${REF_bowtie}.fa
SUFFIX=lang_r2
Sample=${SAMPLE}_${SUFFIX}


#### note that assumes the reads have been trimmed and peared


bowtie2 -p 16 -t --end-to-end --no-unal  -x ${REF_bowtie}  \
	   -U ${OUT_PATH}/${SAMPLE}.trimmed.pear.all.filter.fastq  \
	   -S ${OUT_PATH}/${Sample}.trimmed.pear.all.sam 2> ${OUT_PATH}/${Sample}.bowtie.log 


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
