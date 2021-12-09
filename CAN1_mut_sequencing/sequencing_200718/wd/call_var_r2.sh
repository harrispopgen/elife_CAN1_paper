
FQ_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/raw_reads/nextseq_200718/fastq
OUT_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/sequencing_200718/wd


IFS=$'\n'

for i in `cat sample_seq_info_200718.txt `
do
	sampleseq=`echo $i | cut -f1`
	sample=`echo $i | cut -f6`
	qsub  ../../script/positionFreqs_Lang_new_v2.sh  $sampleseq  ${FQ_PATH}  ${OUT_PATH}  ${sample} 
done




