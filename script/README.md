
- filter vcf files
`process_vcf.sh`

- then call mutations for each individual
- single bp mutations with derived allele freq <0.5 (for fig 1)
`yeast_1002_indi_count_single_main.py`
- rare single bp mutations with with allele count <= 4 (for fig 1)
`yeast_1002_indi_count_single_AC4_main.py`

- single bp mutations with derived allele freq <0.5 with mutations randomly assigned to one haplotype
`yeast_1002_indi_count_single_main_rand.py`

- mutations with triplet context with derived allele freq <0.5
`yeast_1002_indi_count_triplet_main.py`

- synonymous mutations only with derived allele freq <0.5
`yeast_1002_indi_count_single_syn_all_main.py`
