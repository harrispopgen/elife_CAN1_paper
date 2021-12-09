## the pipeline to generate outgroup sites 

## generate all the S. para sites (SNP and non-SNP) which have been aligned to S.cer genome
python yeast_1002_preprocess_sites.py

## generate SNP and indel sites from the alignment
python yeast_1002_preprocess_sites_old.py

## only extract indels 
Rscript extract_indel_filter.R


## merge indels
source merge_indel.sh

## convert chr number in indels to latin
Rscript cvt_chr_indel_bed.R


## exclude indels from mappable sites for each strain
bash SNP_excl_indel.sh

## merge SNP sites from 5 S.para strains (needs large memory: used 8G to run)
Rscript  merge_sites.R
