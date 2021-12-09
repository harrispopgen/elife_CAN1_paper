## calculate for each strain, the C->A ratio 

source("/net/harris/vol1/home/pyjiang/mut_R/plot_PCA_basic.R")



calc_CA_ratio <- function(count_file, outfile)
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
  
  #png("C_A_ratio_dist.png")
  #hist(mut_tb_new_norm[,"C_A"], breaks = 15)
  #dev.off()
  
  C_A_vec <- mut_tb_new_norm[,"C_A"]
  ## simulate a distribution based on C->A ratios in the population
  dist_simu <- sample(C_A_vec, 10000, replace=TRUE)
  
  ## then calculate for each strain the empirical p value 
  p_vec <- NULL
  for(i in 1:nrow(mut_tb_new_norm))
  {
    pi <- sum(dist_simu >C_A_vec[i])/10000
    p_vec <- c(p_vec, pi)
  }
  
  ## assign each individual with its p value and -log(p).
  strains <- rownames(mut_tb_new_norm)
  pval_tb <- data.frame(strains =strains, pval= p_vec )
  pval_tb[,"-logP"] <- -log(pval_tb[,"pval"])
  
  ## write to file 
  write.table(pval_tb, outfile, sep="\t", quote=FALSE, row.names=FALSE)
}

calc_CA_ratio("../data/yeast_single_indi_mut_count_all_AC_4.txt",  "../data/y1002_C_A_ratio_pval.txt")

