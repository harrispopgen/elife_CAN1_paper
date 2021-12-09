## combine all libs that are true single nucleotide variants 

library(data.table)

combine_all_mut <- function()
{
  all_uniq_libs <- read.table("../data/all_libs_uniq_count.txt")
  
  full_mut_tb <- NULL 
  for (lib in all_uniq_libs[,1])
  {
    single_nucl_mut_file <- paste("../all_mut_wd_excl_MNM/", lib, "_final_single_var.txt", sep="")
    single_nucl_mut <- read.table(single_nucl_mut_file, header= TRUE)
    
    single_nucl_mut[,"sample"] <- lib
    single_nucl_mut[,"strain"] <- tstrsplit(single_nucl_mut[,"sample"], "_")[[1]]
    
    full_mut_tb <- rbind(full_mut_tb, single_nucl_mut)
  }
  
  write.table(full_mut_tb, "../all_mut_wd_excl_MNM/all_mut_combined.txt", sep="\t", quote=FALSE, row.names=FALSE)
}

combine_all_mut()
