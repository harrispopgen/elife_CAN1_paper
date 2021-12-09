## This repository contains scripts/pipelines for [Jiang et al 2021. eLife](https://elifesciences.org/articles/68285): A modified fluctuation assay reveals a natural mutator phenotype that drives mutation spectrum variation within Saccharomyces cerevisiae

## fig 1: analyzing mutation spectra from polymorphism data: preprocess of aligning S.paradoxus to S.cerevisiae: Spar_preprocess/, scripts for generating plots are in script/

- need to align outgroup (S.para) to S.cer first.
see directory: `Spar_preprocess/`


- script to plot fig1 and related supplementary figures
`script/plot_fig1_and_suppl.R`
`script/plot_fig1_PCA_Scer_AC4.R`

## fig 3: fluctuation assay: fluc/
`plot_mut_rate_fig3.R`

## fig 4: mutation spectra from de novo mutations: CAN1_mut_sequencing/ 

check the readme file in that directory. first call mutations from raw illumina sequencing reads, identify MNMs, and get the final set of single nucleotide mutations

- script to plot fig 4
`CAN1_mut_sequencing/script/fig4_from_10_21_final.R`

## fig 5: C to A ratio enrichment: dipCBS1782/

refer to the README file in that directory. 

## fig 6: mutation rate estimates from strains transformed with different OGG1 plasmids: OGG1/

figure 6:
`OGG1/plot_mut_rate.R`
