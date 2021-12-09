## extract genes (from bed files) that overlap the sites which have ancestral alleles (also in bed format)

gene_bed_file=/net/harris/vol1/project/yeast_1002/data/R64_1_1_all_gene_coord.bed

gene_gff_file=/net/harris/vol1/project/yeast_1002/data/R64_1_1_all_gene_coord.gff

sites_anc_file=../Spar_preprocess/call_consensus/data/Scer_para_consol_SNPs.txt

source loadBedtools

## to include genes if it has any coverage in ancestral site file 
bedtools intersect  -wa  -a ${gene_bed_file}  -b ${sites_anc_file}  > ../data/gene_coord_in_anc.bed

### will use bcftools intersect to get information of both a and b (will need strand information)

bedtools intersect  -wo -a  ${gene_gff_file}  -b ${sites_anc_file}  > ../data/gene_gff_in_anc.bed

### will extract sub columns to contain only the chr, coord, gene name and strand (similar to Scer_para_consol_SNPs.txt format, with extra columns)
awk 'BEGIN {OFS="\t"}; { print $10,$11,$12,$13,$7}' ../data/gene_gff_in_anc.bed  > ../data/gene_gff_in_anc_consise.bed

### notice that there are certain cases, where at one SNP location, there are more than one ORF, which will result in two strand information
### will keep the SNP if it has consistent strand informaiton, otherwise remove this SNP.


## to generate a gene list that a gene >0.9 coverage are included

Rscript get_orf_in_ance.R

## extract the fasta file (containing all orfs) to the subset of genes 
## (need to run python in a separate window as R; conflict)
python extract_genes_consider.py

