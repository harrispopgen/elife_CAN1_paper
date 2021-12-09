## fig 5 contains three plots, fig 5A and 5B are from 1011 dataset 
## while fig 5C includes incoporating variants called from a newly seqeunced CBS1782 strain with the existing 1011 dataset

## fig 5A:
## simulate how unlikely to see the C to A ratios given the overall distribution of the strains
`script/calc_CA_ratio_pval.R`
`script/calc_CA_ratio_pval_newsimu.R`
## plot fig 5A
`script/plot_fig5A_CA_ratio_tree_pval.R`

## fig 5B:
`script/plot_fig5B_sub_tree.R`

######
## fig 5C
#####
## call raw variants for CBS1782
`script/call_var.sge`

## combine variants from CBS1782 with 1011 yeast genomes
`script/combine_gvcf.sge`

## count mutations from rare polymorphisms given AC cutoff 
`script/record_mut_AC.sge`

## plot 5C
`script/plot_fig5C_CA_ratio_strains.R`
