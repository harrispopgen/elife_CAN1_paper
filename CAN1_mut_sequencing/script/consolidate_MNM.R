## consolidate MNM
library(data.table)

consolidate_MNM <- function(all_MNM_file)
{
  full_tb <- read.table(all_MNM_file, header =TRUE)
  
  pair_mut_tb_all <- NULL 
  
  full_tb[,"file"] <- as.factor(full_tb[,"file"])
  for (file in levels(full_tb[,"file"]))
  {
    subtb <- full_tb [full_tb[,"file"] == file,] 
    
    ## find records that contain either site in the remaining tables (will merge)
    site1 <- subtb[,"pos1"]
    site2 <- subtb[,"pos2"]
    set <- c(site1, site2)
    
    ## if a site is present in only one pair, it should occur once in the "set"
    ## otherwise, it will occur multiple times 
    ## get a set of only one pair of elements 
    duplicated_pos <-  unique(set[duplicated(set)])
    set_vals <- unique(set)
    unique_pos <- setdiff(set_vals, duplicated_pos)
    ## extract rows that are present in unique_pos
    pair_mut_rows <- site1 %in% unique_pos
    
    pair_mut_tb <- subtb[pair_mut_rows,]
    pair_mut_tb_all <- rbind(pair_mut_tb_all, pair_mut_tb)
    
    complex_mut_tb <- subtb[!pair_mut_rows,]
    
    ## deal with complex mut 
    while (nrow(complex_mut_tb)!=0 )
    {
      complex_site1 <- complex_mut_tb[1, "pos1"]
      complex_site2 <- complex_mut_tb[1, "pos2"]
      
      complex_set <- c(complex_site1,complex_site2 )
      
      remain_tb <- complex_mut_tb[2:nrow(complex_mut_tb), ]
      
      ## find all pos that are linked to the initial pair of sites 
      newset <- NULL
      while(!setequal (newset, complex_set))
      {
        if(is.null(newset))
        {
          newset <- complex_set
        }
        
        idx <- (remain_tb[,"pos1"] %in% complex_set )| (remain_tb[,"pos2"] %in% complex_set )
        complex_set <- newset
        newset <- unique(c(remain_tb[idx,"pos1"], remain_tb[idx,"pos2"]))
      }
      
      idx1 <- complex_mut_tb[, "pos1"] %in% complex_set
      ## remaining sites that are not in idx1 
      remain_set <- setdiff(complex_set, unique(complex_mut_tb[idx1, "pos1"]) )
      idx2 <- complex_mut_tb[, "pos2"] %in% remain_set
      
      variable1_tb <- complex_mut_tb [idx1, c("pos1", "ref1", "variable1")]
      names(variable1_tb) <- c("pos", "ref", "variable")
      variable2_tb <- complex_mut_tb [idx2, c("pos2", "ref2", "variable2")]
      names(variable2_tb) <- c("pos", "ref", "variable")
      
      variable_tb <- rbind(variable1_tb, variable2_tb)
      variable_unique_tb <- variable_tb[!duplicated(variable_tb[,"pos"]),]
        
      ## write the complex mutation into a file 
      complex_out_file <- paste("../MNM/final_complex/", file, sep="")
      write.table(variable_unique_tb, complex_out_file, sep="\t", quote=FALSE, row.names=FALSE)
      
      
      ## get the rest of rows 
      complex_mut_tb <- complex_mut_tb[!(idx1 | idx2), ]
      
    }
  }
  
  ## output double mutants of MNM 
  ## generate an extra field of strain 
  pair_mut_tb_all[, "strain"]<- tstrsplit(pair_mut_tb_all[,"file"], split="_")[[1]]
  
  write.table(pair_mut_tb_all,"../MNM/final_MNM_double.txt", sep="\t", quote=FALSE, row.names=FALSE)
  
  consolidate_complex_MNM()
}

## consolidate complext MNM into one file 
consolidate_complex_MNM <- function()
{
  dir <- "../MNM/final_complex/"
  complex_MNM <- data.frame()
  for (file in list.files(dir))
  {
    file_name<- paste(dir, file, sep="")
    tb <- read.table(file_name, header= TRUE)
    mut_str <- ""
    for (i in 1:nrow(tb))
    {
      strl <- paste(tb[i,], collapse ="_")
      if (i==1)
      {
        mut_str <- strl
      }
      else
      {
        mut_str <- paste(mut_str, strl, sep=";")
      }
    }
    new_line <- c(file, mut_str)
    complex_MNM <- rbind(complex_MNM, new_line)
  } 
  names(complex_MNM) <- c("file", "mutation")
  complex_MNM[, "strain"]<- tstrsplit(complex_MNM[,"file"], split="_")[[1]]
  
  write.table(complex_MNM,"../MNM/final_MNM_complex.txt", sep="\t", quote=FALSE, row.names=FALSE )
}



consolidate_MNM("../MNM/final_MNM_list.txt")
