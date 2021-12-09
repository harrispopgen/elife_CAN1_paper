## get candidate MNM mutations 


## the threshold only applies if two sites have exactly one count,
## the relative frequency is ignored if one of the count is >1. 
## read mutation calling file from a single pool (SNV and indel files are separate)
## then pull out sites that are within 10bp from each other. 
## they should also have similar allele frequencies 
get_candidate_MNM <- function(snv_file, indel_file, threshold=0.09)
{
  snv <- read.table(snv_file, header=TRUE, stringsAsFactors = FALSE,
                    colClasses = c("numeric", "character", "character", "numeric",
                                   "character", "numeric", "numeric", "numeric"))
  indel <- read.table(indel_file, header=TRUE, stringsAsFactors = FALSE,
                      colClasses = c("numeric", "character", "character", "numeric",
                                     "character", "numeric", "numeric", "numeric"))
  indel["len"] <- NULL
  
  full <- rbind(snv, indel)
  ## sort by position 
  full_order <- full[order(full[,"pos"]),]
  ## then loop through mutations, to find a pair that is within 10bp of current site 
  candi_pair <- data.frame(pos= numeric(), ref= character(), variable= character(),
                           value= numeric(), count = numeric(), 
                            pos.1= numeric(), ref.1 = character(),
                           variable.1 = character(), value.1= numeric(), count.1= numeric())
  
 
  n <- nrow(full_order)
  for(i in 1:(n-1))
  {
    mut1 <- full_order[i,]
    pos1 <- mut1["pos"]
    j=i+1
    while(j <= n)
    {
      mut2 <- full_order[j,]
      pos2 <- mut2["pos"]
      
      ## if the closest mutation is within 10bp 
      if (pos2 >(pos1+10) )
      {
          break
      }
      else if (pos2 == pos1) ## if two sites are the same, skip (not possible to be present at the same time)
      {
          j=j+1
          next
      }
      else
      {
        ##if they have similar frequencies as well if they have 1 count at each site 
        if (mut1["count"]==1 & mut2["count"]==1)
        {
          frq1 <- mut1["frq"]
          frq2 <- mut2["frq"]
          if(frq1>= frq2)
          {
            margin <- 1- frq2/frq1
          }
          else
          {
            margin <- 1- frq1/frq2
          }
          if (margin <= threshold)
          {
            ## record current pair into list, will output count as well 
            mut_pair <- c(mut1[c(1,3,5,6, 8)], mut2[c(1,3,5,6,8)])
            mut_pair_df <- data.frame(mut_pair)
            candi_pair<- rbind(candi_pair, mut_pair_df) 
          }
        }
        else 
        {
          ## record current pair into list 
          mut_pair <- c(mut1[c(1,3,5,6, 8)], mut2[c(1,3,5,6,8)])
          mut_pair_df <- data.frame(mut_pair)
          candi_pair<- rbind(candi_pair, mut_pair_df) 
        }
        
        j=j+1
      }
      
    }
  }
  names(candi_pair) <- c("pos1", "ref1", "variable1", "read_count1", "mut_count1", "pos2", "ref2", "variable2", "read_count2", "mut_count2")
  return (candi_pair)
}




get_candidate_MNM_tofile <- function(snv_file, indel_file, outfile)
{
  candi_pair <- get_candidate_MNM(snv_file, indel_file)
  if(nrow(candi_pair)>=1)
  {
    write.table(candi_pair, outfile, quote=FALSE, row.names=FALSE)
    return (1)
  }
  return (0)
}

## will also output a table that has the corresponding information needed to validae MNM 
## examine every pool, and output if any pool has any candidate MNMs 
## use filter20 for consistency 
get_candidate_MNM_tofile_all <- function()
{
  all_libs_uniq <- read.table("../data/all_libs_uniq_count.txt")
  # test
  #all_libs_uniq <- read.table("../data/all_libs_uniq_count_GIL.txt")
  
  ## read in the full sample info table 
  sample_info_all <- read.table ("../data/sample_info_all.txt")
  ## test for GIL
  #sample_info_all <- read.table ("../data/sample_info_all_GIL.txt")
  validate_MNM_info_tb <- NULL 
  for (i in 1:nrow(all_libs_uniq))
  {
    pool_prefix <- all_libs_uniq[i,1]
    snp_file<- paste("../all_mut_wd/", pool_prefix, "_filter20_SNV.txt", sep="")
    indel_file<- paste("../all_mut_wd/", pool_prefix, "_filter20_indel.txt", sep="")
    
    out_file <- paste("../MNM/candidate/", pool_prefix, "_MNM_candidate.txt", sep="")
    val <- get_candidate_MNM_tofile(snp_file, indel_file, out_file)
    if(val==1)
    {
      candidate_file <- paste(pool_prefix, "_MNM_candidate.txt", sep="")
      sample_info <- sample_info_all[sample_info_all[,4] == pool_prefix,]
      seq_info <- paste("sequencing_", strsplit(sample_info[,7], split="_")[[1]][4], sep="")
      sample_seq_name <- sample_info[,1]
      strain <- sample_info[,6]
      record <- c(candidate_file, seq_info, sample_seq_name, strain )
      validate_MNM_info_tb <- rbind(validate_MNM_info_tb, record)
    }
  }
  write.table (validate_MNM_info_tb, "../data/MNM_validate_sample_info.txt", sep="\t", quote=FALSE, 
               row.names=FALSE, col.names = FALSE)
}


get_candidate_MNM_tofile_all()

