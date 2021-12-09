declare -a dirs=("sequencing_190610"  "sequencing_190806"  "sequencing_190915" \
	"sequencing_190916"  "sequencing_200121"  "sequencing_200302" \
	"sequencing_200718"  "sequencing_200721"  "sequencing_200918" \
	"sequencing_201026")

for dir in ${dirs[@]}; do
   cd "../"$dir"/wd/"
   Rscript  call_mut_auto_r2.R
   Rscript  get_all_muts_auto.R
   cd ../../script/
done

