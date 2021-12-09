## this is the pipeline of processing natural polymorphisms to count mutation 

## 1. process 1011 Scer vcf 
bash process_vcf.sh

## count single nucleotide mutation (for fig1)
python yeast_1002_indi_count_single_main.py

## count single nucleotide mutation, for variants that are shared by multiple haplotype, only count once in a random individual
python yeast_1002_indi_count_single_main_rand.py

## count single nucleotide mutation, syn only
## first generate a bed file that include the regions considered (all ORFs, gene has >0.9 coverage)
bash extract_genes_considered.sh
Rscript keep_consistent_ORF_SNP.R
python yeast_1002_indi_count_single_syn_all_main.py

## count single nucleotide mutation, in triplet context
python yeast_1002_indi_count_triplet_main.py

## count single nucleotide mutation, with derived allele count <=4 (rare variants, for fig1)
python yeast_1002_indi_count_single_AC4_main.py


## for fig1 suppl figs
Rscript plot_fig1_and_suppl.R

## fig1 rare variant PCA:
Rscript plot_fig1_PCA_Scer_AC4.R
