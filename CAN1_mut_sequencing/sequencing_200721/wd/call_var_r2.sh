
FQ_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/raw_reads/nextseq_200302
OUT_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/sequencing_200302/wd
#REF=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/S288C_CAN1_Lang_primer.fasta


IFS=$'\n'

for i in `cat sample_seq_info_200302.txt `
#for i in `cat GIL_seq_info_200302.txt `
do
	sampleseq=`echo $i | cut -f1`
	sample=`echo $i | cut -f6`
	qsub  ../../script/positionFreqs_Lang_new_v2.sh  $sampleseq  ${FQ_PATH}  ${OUT_PATH}  ${sample} 
done




