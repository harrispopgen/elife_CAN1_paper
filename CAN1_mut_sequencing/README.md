#####
call mutations from each pooled sequencing
#####

each directory beginning with "sequencing_*/wd" contains pipeline to mutations from raw reads from each sequencing batch
run `call_var_S288ref.sh` in each directory, to align reads to the CAN1 in the reference genome first. 

run `script/get_fix_sites_all_samples.sh` to get all fixed sites from each sample, to check with SNPs known in CAN1 for each strain.

then go to each "sequencing_*/wd" directory, to run `call_var_r2.sh`, to align reads to each CAN1 sequence from its strain. and generate a frequency table of each variant at each site. 

run `script/call_mut_in_all_samples.sh` to call mutations. after that, raw mutations are output in the "all_mut_wd" directory. 

#####
analyze MNMs, first find candidates, then combine by looking at aligned reads, then combine complex MNMs if they overlap variants in the same pool  
#####

run `script/MNN_pipeline.sh`


####
excluding variants in MNM, and fitting remaining reads to regenerate counts for SNVs, combine all libs after excluding MNM
###
run `script/refine_SNV.R`
run `script/combine_all_libs_excl_MNM.R`


####
script to generate fig 4 and relavant supplementary figures
####

`script/fig4_from_10_21_final.R`

####
script to calculate numbers of each variant types in the text
####
`script/analyze_mut_type_CAN1.R`
