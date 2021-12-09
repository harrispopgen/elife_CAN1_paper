## usage: ./align_to_Scer.sh  species_shortname  species_fasta

## use lastz to align scaffolds from other genomes to Scer

## give the parameter with the query's fasta file name 
## note: the outgroup S.para have the same chr numbers as S.cer

outgroup_name=$1 ## for example: "CBS432"

S_cere_prefix="/net/harris/vol1/project/yeast_1002/nobackup/S_cere_ref/S288C_R64-1-1_chr/chr"

S_other_prefix="/net/harris/vol1/project/yeast_1002/nobackup/S_para/"

## declare an array variable
declare -a chrs=("I" "II" "III" "IV" "V" "VI" "VII" "VIII" "IX" "X" "XI" "XII" "XIII" "XIV" "XV" "XVI")

for chr in "${chrs[@]}"
do
	echo "align chr$chr"
	S_cere_chr="$S_cere_prefix$chr.fa"
	S_other="$S_other_prefix/${outgroup_name}.genome.chr$chr.fa"
	lastz  $S_cere_chr $S_other  --format=lav  --traceback=500M   --output="../data/S_cer_${outgroup_name}_chr$chr.lav"
done


