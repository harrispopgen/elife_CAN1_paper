## this script align to S288 reference (will double check with SNP first, before calling more accurate variants)

## path for the fastq
FQ_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/raw_reads/nextseq_190610
## output directory (this folder)
OUT_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/sequencing_190610/wd
## reference seq of CAN1_Lang fragment
REF=../../data/S288C_CAN1_Lang_primer.fasta
suffix=001

IFS=$'\n'

for i in `cat sample_seq_info_190610.txt `
do
	sampleseq=`echo $i | cut -f1`
	qsub  ../../script/positionFreqs_Lang_new_gz_v2.sh  $sampleseq  ${FQ_PATH}  ${OUT_PATH}  ${REF} ${suffix}
done




