## process the original vcf to the input for python

## under centos7

module load tabix/0.2.6
## load bcftolls
module load htslib/1.9
module load bcftools/1.9
## load VCFtools
module load VCFtools/0.1.14
## load bedtools
module load bedtools/2.27.1


## the original vcf from 1011 dataset 
vcf=/net/harris/vol1/project/yeast_1002/nobackup/1011Matrix.vcf.gz
outdir=../data/vcf
prefix=1011_SNP_bi_mm0.8_GT


## set ID as chr:pos and get biallelic SNPs
## chr names are 1,2,3...
bcftools annotate --output-type v --remove ID  --set-id +'%CHROM:%POS' $vcf \
	 | bcftools view -m2 -M2  -v snps  | \
     sed 's/chromosome//g' | bgzip -c > $outdir/1011Matrix_SNPs_only_ID_chr_no_biallele.recode.vcf.gz


## filter out genotypes  with missing values > 0.2
## keep only genotype
vcftools --gzvcf $outdir/1011Matrix_SNPs_only_ID_chr_no_biallele.recode.vcf.gz \
    --max-missing 0.8  --recode  --recode-INFO-all --stdout \
    | bcftools  annotate -x 'FORMAT'  \
        -o $outdir/1011Matrix_SNPs_only_ID_chr_no_biallele_mm0.8.GT.vcf


######
## note: needs to chang the reference of the outgroup so that it filter out the indel and merge SNPs from Spara
######
## rename chr from numbers to latin
python -m replace_chromosome_names  --cols 1  --comment-char "#"  -m ../data/yeast_chr_no_to_latin.txt \
	-o $outdir/${prefix}.vcf \
	$outdir/1011Matrix_SNPs_only_ID_chr_no_biallele_mm0.8.GT.vcf


bgzip   $outdir/${prefix}.vcf

tabix -p vcf  $outdir/${prefix}.vcf.gz



## 2. keep SNPs in uniquely mapped regions
## and exclude repeat masked regions

bedtools  intersect  -u  -a  $outdir/${prefix}.vcf.gz \
	 -b ../data/uregions_100_2011_tb.bed  \
         -header | bedtools intersect -v  -a stdin \
	 -b  ../data/repeat_etc_mask.gff  \
         -header | bgzip -c >  ../data/vcf/${prefix}_u_rm.vcf.gz

tabix -p vcf  ../data/vcf/${prefix}_u_rm.vcf.gz




###########
## use all individuals excluding close relatives for AC count
## Note!!: need to update AC, AN count !!
###########
bcftools view -S  /net/harris/vol1/project/mut_survey/species/Scer/beer/data/1011_strains_no_close.txt --min-ac=1 \
         ../data/vcf/${prefix}_u_rm.vcf.gz  \
		 | bcftools  +fill-tags  -- -t AN,AC,AF \
        | bgzip -c  > ../data/vcf/${prefix}_u_rm_final_no_close.vcf.gz



## split by chr

declare -a chrs=( "chrI" "chrII" "chrIII" "chrIV" "chrV" "chrVI" "chrVII" "chrVIII" "chrIX" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI" )
chr_counter=1
for i in "${chrs[@]}"; do vcftools  --gzvcf   ../data/vcf/${prefix}_u_rm_final_no_close.vcf.gz \
         --chr $i  --recode --recode-INFO-all \
         --stdout | bgzip -c >  ../data/vcf/${prefix}_u_rm_final_no_close.chr${chr_counter}.vcf.gz ;
	chr_counter=$((chr_counter+1))
done



#########
## generate annotated vcf for analysis of only syn mutations
#########


## the annotated vcf is generated from /net/harris/vol1/home/pyjiang/software/snpEff/yeast_1002_snpEff2.sge

bgzip  $outdir/${prefix}.annot.vcf

tabix -p vcf  $outdir/${prefix}.annot.vcf.gz

#### follow the same procedure as the unannotated

## 2. keep SNPs in uniquely mapped regions
## and exclude repeat masked regions

bedtools  intersect  -u  -a  $outdir/${prefix}.annot.vcf.gz \
         -b ../data/uregions_100_2011_tb.bed  \
         -header | bedtools intersect -v  -a stdin \
         -b  ../data/repeat_etc_mask.gff  \
         -header | bgzip -c >  ../data/vcf/${prefix}_u_rm.annot.vcf.gz


tabix -p vcf  ../data/vcf/${prefix}_u_rm.annot.vcf.gz


###########
## use all individuals excluding close relatives for AC count
###########
bcftools view -S  ../data/1011_strains_no_close.txt --min-ac=1 \
         ../data/vcf/${prefix}_u_rm.annot.vcf.gz  \
		 | bcftools  +fill-tags  -- -t AN,AC,AF \
        | bgzip -c  > ../data/vcf/${prefix}_u_rm_final_no_close.annot.vcf.gz


## split by chr

declare -a chrs=( "chrI" "chrII" "chrIII" "chrIV" "chrV" "chrVI" "chrVII" "chrVIII" "chrIX" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI" )
chr_counter=1
for i in "${chrs[@]}"; do vcftools  --gzvcf   ../data/vcf/${prefix}_u_rm_final_no_close.annot.vcf.gz \
         --chr $i  --recode --recode-INFO-all \
         --stdout | bgzip -c >  ../data/vcf/${prefix}_u_rm_final_no_close.annot.chr${chr_counter}.vcf.gz ;
        chr_counter=$((chr_counter+1))
done
