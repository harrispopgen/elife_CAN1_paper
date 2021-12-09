## validate MNM candidates 

IFS=$'\n'

for i in `cat ../data/MNM_validate_sample_info.txt `
do
  	qsub validate_MNM_single.sge  $i 20
done


for i in `cat ../data/MNM_validate_sample_info.txt `
do
	qsub validate_MNM_single.sge  $i 40
done 

