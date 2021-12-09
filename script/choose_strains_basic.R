## write a function to randomly choose strains 


## annot_tb_new: tb with annotation
## num_per_group:  maximum num of strain to subsample per group 
## strain_list: strain list to choose from (before subsampling)
subsample_strains <- function(annot_tb, num_per_group, strain_list)
{
  ## only include strains that are in strain_list (if not NULL)
  if (!is.null(strain_list))
  {
    annot_tb_new <- annot_tb [ annot_tb[,"ID"] %in% strain_list, ]
  }
  
  new_tb <- NULL
  
  annot_tb_new[,"group"] <- factor(annot_tb_new[,"group"] )
  
  for (group in levels(annot_tb_new[,"group"]) )
  {
    subtb <- annot_tb_new [annot_tb_new[,"group"] == group,]
    if (nrow(subtb) <= num_per_group)
    {
      new_tb <- rbind(new_tb, subtb)
    }
    else
    {
      ## randomly choose 30, without replacement
      sample_ids <- sample.int(nrow(subtb), num_per_group)
      new_tb <- rbind(new_tb, subtb[sample_ids,])
    }
  }
  return (new_tb)
}



## this is the original function 
## strain_list is speficied if a subset of the original strain is includeed 
choose_strains_rand <- function(annot_file, strains_exclude_file, num_per_group, strain_list)
{
  annot_tb <- read.table(annot_file, sep=",", header=TRUE, stringsAsFactors = FALSE)
  
  ## first exclude the introgressed strains 
  excl_strains <- read.table(strains_exclude_file, stringsAsFactors = FALSE)
  
  
  ## after remove introgressed strains 
  annot_tb_new <- annot_tb [! annot_tb[,"ID"] %in% excl_strains[,1], ]
  
  new_tb <- subsample_strains(annot_tb_new, num_per_group, strain_list)
  
  return (new_tb)
}


## exclude strains from certain groups
choose_strains_rand2 <-  function(annot_file, strains_groups_excl, num_per_group, strain_list)
{
  
  annot_tb <- read.table(annot_file, sep=",", header=TRUE, stringsAsFactors = FALSE)
  
  ## after remove introgressed strains 
  annot_tb_new <- annot_tb [! annot_tb[,"group"] %in% strains_groups_excl, ]
  
  new_tb <- subsample_strains(annot_tb_new, num_per_group, strain_list)
  
  return (new_tb)
  
}
