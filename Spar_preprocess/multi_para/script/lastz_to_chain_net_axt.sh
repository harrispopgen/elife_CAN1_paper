## note that since maf to axt needs more information (each fasta record name)
## would switch to lav to psl format instead

## convert the lav format to axt format

## 
echo "lav to psl"
outgroup_name=$1


wd_dir="../data"

S_other_prefix="/net/harris/vol1/project/yeast_1002/nobackup/S_para/"

S_other_fa="${S_other_prefix}/${outgroup_name}.genome.fa"

Scer_ref_dir="/net/harris/vol1/data/Pengyao/nobackup/S288C_ref"

declare -a chrs=("I" "II" "III" "IV" "V" "VI" "VII" "VIII" "IX" "X" "XI" "XII" "XIII" "XIV" "XV" "XVI")

for chr in "${chrs[@]}"
do
    lav_name="$wd_dir/S_cer_${outgroup_name}_chr$chr.lav"
    lavToPsl  $lav_name    "$wd_dir/S_cer_${outgroup_name}_chr$chr.psl"
done

## generate 2bit file for scaffold

echo "fa to 2bit"
faToTwoBit  $S_other_fa  "$wd_dir/${outgroup_name}.2bit"

## psl to chain
echo "to chain"
for chr in "${chrs[@]}"
do
	axtChain -minScore=10000 -linearGap=medium  "$wd_dir/S_cer_${outgroup_name}_chr$chr.psl" \
		 "$Scer_ref_dir/sacCer3.2bit" \
		 "$wd_dir/${outgroup_name}.2bit"   -psl  stdout |
	chainAntiRepeat "$Scer_ref_dir/sacCer3.2bit"  \
		 "$wd_dir/${outgroup_name}.2bit"   stdin   "$wd_dir/S_cer_${outgroup_name}_chr$chr.chain"
done


## merge chain and sort
echo "merge chain"
ls $wd_dir/S_cer*${outgroup_name}*chain | chainMergeSort  -inputList=stdin | gzip -c > "$wd_dir/sacCer3.${outgroup_name}.all.chain.gz"

## calculate the size of S_other genome (per chr)
faSize  $S_other_fa  -detailed > "$wd_dir/${outgroup_name}.sizes"

chainPreNet "$wd_dir/sacCer3.${outgroup_name}.all.chain.gz" \
	  "$Scer_ref_dir/sacCer3.chrom.sizes"  \
	  "$wd_dir/${outgroup_name}.sizes"  "$wd_dir/sacCer3.${outgroup_name}.all.prenet.chain"


## net chains
echo "net chains"
chainNet "$wd_dir/sacCer3.${outgroup_name}.all.prenet.chain" -minSpace=1 \
	  "$Scer_ref_dir/sacCer3.chrom.sizes" \
	  "$wd_dir/${outgroup_name}.sizes"  stdout /dev/null \
	  | netSyntenic stdin  "$wd_dir/noClass_${outgroup_name}.net"

## net to axt (note: needs to use axtSort, otherwise, the axt are not ordered from small to large in coordinate)
echo "net to axt"
netToAxt  "$wd_dir/noClass_${outgroup_name}.net" \
		 "$wd_dir/sacCer3.${outgroup_name}.all.prenet.chain" \
		 "$Scer_ref_dir/sacCer3.2bit" \
		 "$wd_dir/${outgroup_name}.2bit"   stdout | \
		  axtSort stdin  "$wd_dir/Scer_${outgroup_name}.sorted.net.axt"


## convert axt to maf (to load at UCSC browser)
axtToMaf "$wd_dir/Scer_${outgroup_name}.sorted.net.axt" \
	"$Scer_ref_dir/sacCer3.chrom.sizes" \
	"$wd_dir/${outgroup_name}.sizes" \
	"$wd_dir/Scer-${outgroup_name}.maf" -tPrefix=sacCer3.

