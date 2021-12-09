## by excluding variants in MNM 
## and fitting remaining reads to regenerate counts 


source("/net/harris/vol1/home/pyjiang/mut_R/count_mut_pool.R")




refine_SNV <- function(sample_prefix_file, mut_dir, outdir)
{
  sample_prefix <- read.table(sample_prefix_file)
  
  MNM_sites_dir <-"../MNM/final_MNM_exclude"
  all_MNM_sites_files <- list.files(MNM_sites_dir)
  
  MNM_remain_single_dir <- "../MNM/final_single"
  MNM_remain_single_files <- list.files(MNM_remain_single_dir)
  
  for(prefix in sample_prefix[,1])
  {
    ## read in SNV file and indel file separately
    SNV_file <- paste(mut_dir, prefix, "_filter40_SNV.txt", sep="")
    indel_file <- paste(mut_dir, prefix, "_filter20_indel.txt", sep="")
    
    SNV_tb <- read.table(SNV_file, header=TRUE)
    SNV_tb[,"len"] <- 1
    indel_tb <- read.table(indel_file, header=TRUE)
    
    full_tb <- rbind(SNV_tb, indel_tb)
    full_tb[,"var"] <- paste(full_tb[,"pos"], full_tb[,"variable"],sep="_")
    
    ## then remove sites if they are present in MNM
    MNM_site_file_name <- paste(prefix, "_MNM_candidate.txt", sep="" )
    ## exclude MNM sites if have any 
    if (MNM_site_file_name %in% all_MNM_sites_files)
    {
      MNM_site_file_name_full <- paste("../MNM/final_MNM_exclude/", MNM_site_file_name, sep="")
      MNM_sites <- read.table(MNM_site_file_name_full, header =TRUE)
      
      MNM_sites[,"var"] <- paste(MNM_sites[,"pos"], MNM_sites[,"variable"],sep="_")
      
      idx <- full_tb[,"var"] %in% MNM_sites[,"var"]
      full_tb <- full_tb[!idx,]
    }
    
    ## for the sites that overlap with an MNM site, change the remaining read counts 
    if (MNM_site_file_name %in% MNM_remain_single_files)
    {
      single_remain_file <- paste("../MNM/final_single/", MNM_site_file_name, sep="")
      single_sites <- read.table(single_remain_file, header =TRUE,
                                 colClasses = c("numeric",  "character", "character", "numeric"))
      
      single_sites[,"var"] <- paste(single_sites[,"pos"], single_sites[,"variable"],sep="_")
      
      idx <- full_tb[,"var"] %in% single_sites[,"var"]
      full_tb_rest <- full_tb[!idx,]
      full_tb_adjust <- full_tb[idx,]
      ## change the read count into the ones 
      for(i in 1:nrow(full_tb_adjust))
      {
        var <- full_tb_adjust[i,"var"]
        full_tb_adjust[i,"value"] <- single_sites[single_sites[,"var"] == var,"read_count"]
        full_tb_adjust[i,"frq"] <- full_tb_adjust[i,"value"] / full_tb_adjust[i,"reads_all"]
      }
      full_tb <- rbind(full_tb_rest, full_tb_adjust)
    }
    
    full_tb[,"var"] <- NULL
    
    ## re-estimate the mutant count 
    full_tb <- count_mut_tb_core(full_tb)
    
    ##############
    ## keep all indels now 
    ##############
    ## only include variants of sinlge bp 
    #full_tb <- full_tb [full_tb[,"len"] == 1,] 
    
    outfile <- paste(outdir, prefix, "_final_single_var.txt", sep="")
    write.table(full_tb, outfile, quote=FALSE, row.names=FALSE,sep="\t")
  }
}

refine_SNV("../data/all_libs_uniq_count.txt", "../all_mut_wd/", "../all_mut_wd_excl_MNM/")
