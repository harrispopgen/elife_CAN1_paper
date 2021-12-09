## script to get the final list of strains used in the paper 


source("/net/harris/vol1/home/pyjiang/mut_R/choose_strains_basic.R")


## read in the strains that already exclude close relatives
## then exclude strains that have extensive introgressions from S.para (clade 2,9,10)
choose_strains <- function()
{
  #### note that i introduced the randomness of PCA by randomly choosing 30 strains every time 
  #### which 
  strains <- read.table("/net/harris/vol1/project/mut_survey/species/Scer/beer/data/1011_strains_no_close.txt")
  
  
  new_tb <- choose_strains_rand2("/net/harris/vol1/project/yeast_1002/data/strain_anno_group.txt",
                                c(2,9,10), 30, strains[,1])
  
  ## write both annotation file and strain name only file
  write.table(new_tb,"../data/strain_subset_final_AC4_anno_no_close.csv", sep=",", quote=FALSE, row.names=FALSE)
  write.table(new_tb[,"ID"],"../data/strain_subset_final_AC4_no_close.txt", quote=FALSE, row.names=FALSE)
}


## only choose strains if a file does not exist 
if(!file.exists("../data/strain_subset_final_AC4_no_close.txt"))
{
  choose_strains()
}




