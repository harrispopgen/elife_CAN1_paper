## calculate p value from new simulations (randomly sample
## mutations based on the same number of AC for each strain)

## calculate for each strain, the C->A ratio 

source("/net/harris/vol1/home/pyjiang/mut_R/plot_PCA_basic.R")

calc_CA_ratio_pval <- function(count_file, simu_dir_prefix, outfile)
{
  ## read in individual mutation count for each strain with AC cutoff=4
  all_count <- read.table(count_file, header =TRUE)
  all_count_mut <- all_count[-c(13,14), ] ## exclude indels 
  
  ## exclude strains  which have total mutation count <8
  
  count_tb <- all_count_mut[,2:ncol(all_count_mut)]
  total <- colSums(count_tb)
  all_count_no_outlier <- count_tb[,total >=8 ]
  outlier <- names(count_tb[,total <8 ])
  all_count_no_outlier_tb <- cbind("mut"=all_count_mut[,1] ,all_count_no_outlier )
  
  
  mut_tb_new_norm <- PCA_preprocess(all_count_no_outlier_tb,0, 1,0)
  
  ## calculate CA ratio observed in each strain 
  C_A_vec <- mut_tb_new_norm[,"C_A"]
  strains <- names(C_A_vec)
  pval_vec <- c()
  
  for (i in 1:length(C_A_vec))
  {
    strain <- strains[i]
    ## then read the simulated C_A ratios for each strain 
    simu_CA_ratio_file <- paste(simu_dir_prefix, strain, "_CA_ratio.txt", sep="")
    simu_CA_ratio <- read.table(simu_CA_ratio_file, header=TRUE)
    
    pval <- sum(simu_CA_ratio >= C_A_vec[i])/10000
    pval_vec <- c(pval_vec, pval)
  }
  
  ## output strains with their empircal p-val
  pval_tb <- data.frame(strains =strains, pval= pval_vec )
  pval_tb[,"-logP"] <- -log(pval_tb[,"pval"])
  
  ## write to file 
  write.table(pval_tb, outfile, sep="\t", quote=FALSE, row.names=FALSE)
}

calc_CA_ratio_pval("../../data/Scer_single_no_close_AC4_indi_mut_count_strain.txt", "../data/simu/", "../data/y1002_C_A_ratio_pval_newsimu.txt")
