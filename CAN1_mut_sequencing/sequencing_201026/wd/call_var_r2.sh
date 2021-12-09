
FQ_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/raw_reads/nextseq_201026/fastq
OUT_PATH=/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/sequencing_201026/wd


IFS=$'\n'

for i in `cat sample_seq_info_201026.txt `
do
	sampleseq=`echo $i | cut -f1`
	sample=`echo $i | cut -f6`
	qsub  ../../script/positionFreqs_Lang_new_v2.sh  $sampleseq  ${FQ_PATH}  ${OUT_PATH}  ${sample} 
done




