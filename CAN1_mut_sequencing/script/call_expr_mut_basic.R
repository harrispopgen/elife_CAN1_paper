library(reshape2)


library(data.table)

library(seqinr)

#library(Ckmeans.1d.dp)

source("/net/harris/vol1/home/pyjiang/mut_R/mut_basic_funcs.R")


## given the per site mutation count table, will output which variant type
## 10.11 change ratio to threshold
get_mut_type <- function(can1_mut_candidate, threshold)
{
  ## convert to long table first
  can1_mut_candidate_long <- reshape2::melt(can1_mut_candidate, id.vars=c("pos","chrom","ref","reads_all"))

  ## exclude the sites which are the same as ref
  can1_mut_candidate_long[,"variable"] <- as.character(can1_mut_candidate_long[,"variable"])
  can1_mut_candidate_long_new <- can1_mut_candidate_long [ can1_mut_candidate_long[,"variable"] != can1_mut_candidate_long[,"ref"] ,]

  ## then caclulate frequency of each mutation type (del, ins, A, C, T, G)
  can1_mut_candidate_long_new[,"frq"]<- can1_mut_candidate_long_new[,"value"]/can1_mut_candidate_long_new[,"reads_all"]

  can1_mut_candidate_long_final <- can1_mut_candidate_long_new[can1_mut_candidate_long_new[,"frq"] > threshold, ]

  return(can1_mut_candidate_long_final)
}


## modify to include indels
## function to screen for candidate mutants
## add begin and end parameter (in case the reference is larger than the amplicon)
## begin and end are both 1-based coord, inclusive.
## modify on 10.11, change ratio to 0.65*exp_mut_freq or 0.02, whichever is larger
get_candicate_mut <- function(can1_site_seq_tb, exp_mut_freq, begin=-1,end=-1)
{
  
      ## set filters
      ## 1. only extract sites where mismatch counts > (0.65 * exp_mut_freq)
      
      ratio =0.65
      threshold <- max(ratio * exp_mut_freq, 0.02)
      
      filter_frq_cutoff_SNP <- (can1_site_seq_tb[,"mismatches"]/can1_site_seq_tb[,"reads_all"] ) > threshold
      
      ## include indels
      filter_frq_cutoff_del <- (can1_site_seq_tb[,"deletions"]/can1_site_seq_tb[,"reads_all"] ) > threshold
      filter_frq_cutoff_ins <- (can1_site_seq_tb[,"insertions"]/can1_site_seq_tb[,"reads_all"] ) > threshold
      
      ## 2. exclude sites which are fixed in this line (but differ from ref)
      filter_excl_anc <- (can1_site_seq_tb[,"mismatches"]/can1_site_seq_tb[,"reads_all"] ) < 0.95
      
      ## 3. exclude sites with coverage < 200
      filter_cover <- can1_site_seq_tb[,"reads_all"] >= 200
      
      ## 4. exclude sites that have coverage < 40% of avg coverage at this locus
      mean_cov <- mean(can1_site_seq_tb[,"reads_all"])
      filter_cover2 <- can1_site_seq_tb[,"reads_all"] >= 0.4 * mean_cov
      
      
      can1_mut_candidate <- can1_site_seq_tb[ (filter_frq_cutoff_SNP |filter_frq_cutoff_del |  filter_frq_cutoff_ins) & filter_excl_anc & filter_cover & filter_cover2,
                                              c("pos", "chrom", "ref", "reads_all", "deletions","insertions", "A","C","T","G")]
      
      can1_mut_candidate_long_final <- get_mut_type(can1_mut_candidate, threshold)
      
      if (begin!=-1 | end!=-1)
      {
        can1_mut_candidate_long_final <- can1_mut_candidate_long_final[(can1_mut_candidate_long_final[,"pos"] >= begin) & (can1_mut_candidate_long_final[,"pos"] <= end),]
      }
    
    
    return (can1_mut_candidate_long_final)
}


## get fixed mut sites (alt allele frq > 0.95)
get_fixed_diff_sites <- function(can1_site_seq_tb)
{
  fixed_sites <- can1_site_seq_tb [((can1_site_seq_tb[,"mismatches"]/can1_site_seq_tb[,"reads_all"] ) > 0.95)  & (can1_site_seq_tb[,"reads_all"] >= 200), c("pos","chrom","ref","reads_all", "deletions","insertions", "A","C","T","G")]

  fixed_sites_long_final <- get_mut_type(fixed_sites,0.95)
  
  ## only report fixed sites in 109-1879 (CAN1 coding region)
  fixed_sites_long_final_cod <- fixed_sites_long_final[(fixed_sites_long_final[,"pos"]>=109) & (fixed_sites_long_final[,"pos"]<=1879) ,]

  ## sort it by pos
  fixed_sites_long_final_sort <- fixed_sites_long_final_cod[order(fixed_sites_long_final_cod[,"pos"]),]

  return (fixed_sites_long_final_sort)
}

## replace fixed sites with the fixed sites of this specific strain
get_fixed_diff_sites_tofile <- function(mut_tb_prefix, ref_file, strain_name, out_fa_dir)
{
  mut_tb_file <-  paste(mut_tb_prefix, "lang.filter.all.stats.txt", sep="_")
  tb1<- read.table(mut_tb_file,header=TRUE)
  fixed_sites <- get_fixed_diff_sites(tb1)

  ref <- toupper(read.fasta(ref_file)[[1]])
  for (i in 1:nrow(fixed_sites))
  {
    pos <- fixed_sites[i,"pos"]
    base <- fixed_sites[i, "variable"]
    ref[pos] <- base
  }

  header <- paste(strain_name, "CAN1_Lang", sep="_")

  out_fa_name <- paste(out_fa_dir, strain_name, "_CAN1_Lang.fa", sep="")
  write.fasta(list(ref), list(header),file.out =  out_fa_name)
}

## will get fixed sites for each sample (to check the correctness) of that strain
## to make sure it did not get mixed up
## will output into one table with two columns, col 1 is the library name, col 2 is a string with all fixed mut, separated by ","
get_fixed_diff_sites_for_samples <- function(sample_info_prefix)
{
  sample_info_tb <- paste(sample_info_prefix, ".txt" ,sep="")
  ## read in the sample info  table
  sample_info <-  read.table(sample_info_tb, stringsAsFactors = FALSE)

  nsamples <- nrow(sample_info)

  fix_info_tb <- NULL
  for(i in 1:nsamples)
  {
    input_file <- paste(sample_info[i,1], "lang.filter20.all.stats.txt", sep="_")

    tb <- read.table(input_file, header = TRUE)

    fixed_site_tb <- get_fixed_diff_sites(tb)

    fixed_site_tb[,"fix"] <- paste(fixed_site_tb[,"pos"], fixed_site_tb[,"variable"], sep="")

    fixed_site_str <- paste(fixed_site_tb[,"fix"] , collapse=",")

    fix_info <- c(sample_info[i,1], fixed_site_str)
    fix_info_tb <- rbind(fix_info_tb, fix_info )
  }
  ## write fixed mut to file
  outfile <- paste( sample_info_prefix, "_fix_sites.txt", sep="")
  write.table(fix_info_tb, outfile, quote =FALSE, row.names =  FALSE, col.names =  FALSE, sep="\t")
}



## did not change the coordinate, but just include the variants in the coding region
## the table will never be empty
call_mut <- function(ref_begin, amplicon_begin, amplicon_end, exp_mut_freq, mut_tb_file, outfile)
{
	begin = amplicon_begin -ref_begin+1
	end = amplicon_end -ref_begin+1


	tb1<- read.table(mut_tb_file,header=TRUE)
	#print(floor(mean(tb1[,"reads_all"])))

	can1_candi_mut1<- get_candicate_mut(tb1, exp_mut_freq, begin, end)

	can1_mut_sort <- can1_candi_mut1[order(as.numeric(as.character(can1_candi_mut1[,"pos"]))),]

	write.table(can1_mut_sort, outfile,quote=FALSE,row.names=FALSE,sep="\t")

	## will return the mean coverage of this gene
	return (floor(mean(tb1[,"reads_all"])))


}

## modified on 1.30.20, add prefix r2
## call mut from a table with info for different sample libraries
## sample_info_tb is a txt file which has the following columns:
## libarary_name(after sequence), #of mutants in a pool, coef to normalize (1.5 for O.N., other 1), libaray name without S1 etc.
## modified on 10.8, call mut for both filter20 and filter40
call_mut_all <- function(sample_info_prefix, outdir= "../../all_mut_wd/")
{

  sample_info_tb <- paste(sample_info_prefix, ".txt" ,sep="")
  ## read in the sample info  table
  sample_info <-  read.table(sample_info_tb)

  nsamples <- nrow(sample_info)

  ref_begin <- 31588
  amplicon_begin <- 31694
  amplicon_end <- 33466

  coverage_list <- NULL
  for(i in 1:nsamples)
  {
    exp_mut_freq <- 1/sample_info[i,2]/sample_info[i,3]

    ## for filter MQ 20
    input_file <- paste(sample_info[i,1], "lang_r2.filter20.all.stats.txt", sep="_")

    output_file <- paste(outdir, sample_info[i,4], "_filter20.txt", sep="")

    coverage <- call_mut(ref_begin, amplicon_begin, amplicon_end, exp_mut_freq, input_file, output_file)

    coverage_list <- c(coverage_list, coverage)


    ## for filer MQ 40
    input_file <- paste(sample_info[i,1], "lang_r2.filter40.all.stats.txt", sep="_")

    output_file <- paste(outdir, sample_info[i,4], "_filter40.txt", sep="")

    call_mut(ref_begin, amplicon_begin, amplicon_end, exp_mut_freq, input_file, output_file)

  }

  ## write coverage to table
  cov_tb <- data.frame(sample =sample_info[,4], coverage = coverage_list)
  cov_out <- paste(sample_info_prefix, "_cov.txt" ,sep="")
  write.table(cov_tb,cov_out, sep=" ", quote=FALSE, row.names= FALSE, col.names=FALSE )
}


## count the sites did not pass the coverage criteria
## add new arguments (only considering sites within certain region)
sites_miss_cover <- function (can1_site_seq_file1, can1_site_seq_file2, begin=-1, end=-1)
{
  can1_site_seq_tb1 <- read.table(can1_site_seq_file1, header=TRUE)
  can1_site_seq_tb2 <- read.table(can1_site_seq_file2, header=TRUE)
	sites <- can1_site_seq_tb1 [(can1_site_seq_tb1[,"reads_all"] <= 200) | (can1_site_seq_tb2[,"reads_all"] <= 200), ]

	if (begin!=-1 | end!=-1)
	{
	  sites <-  sites[(sites[,"pos"] >= begin) & (sites[,"pos"] <= end),]
	}
	return (sites)
}


## updated on 3.9.20, add dir
## combine mutations from list of files
## now separated SNV from indel, will pass a parameter as suffix
combine_muts <- function(dir, prefix_list, suffix)
{
  full_mut_tb <- NULL
  for(prefix in prefix_list)
  {
    pp <- paste(dir, prefix, sep="")
    mut_file_name <- paste(pp, suffix, sep="_")
    mut_tb <- read.table(mut_file_name, header=TRUE, sep="\t", stringsAsFactors=FALSE, colClasses="character" )
    if (nrow(mut_tb) >0)
    {
      mut_tb[,"sample"] <- prefix

      full_mut_tb <- rbind(full_mut_tb, mut_tb)
    }
  }
  return (full_mut_tb)
}

## combine mutations from previous mutation file with the newly generated mutation table
## note: previous file has the strain field
## write the combined table into output file
combine_muts_w_file <- function(prev_mut_file, current_mut_tb, outfile)
{
  prev_mut <- read.table(prev_mut_file, header=TRUE)
  all_mut <- rbind (prev_mut, current_mut_tb)

  write.table(all_mut, outfile, sep="\t", quote=FALSE, row.names=FALSE)
}


## get triplet context for each mutant
## here the ref is an amplicon
## suppose the fasta seq is the first
get_context <- function (ref_file, mut_tb_file)
{
  ref <- toupper(read.fasta(ref_file)[[1]])
  mut_tb <- read.table(mut_tb_file, header=TRUE)

  ## for each mut, get its 5' and 3' seq
  for ( i in 1:nrow(mut_tb))
  {
    pos <- mut_tb[i,"pos"]
    ## 5' and 3' seq of the mut
    context <- paste(ref[pos-1], ref[pos+1], sep="")
    mut_tb[i,"context"] <- context
  }
  return (mut_tb)
}




## modified on 7.8.19, to include insertion and deletion counts
## this function is only suitable to count mutation spec that are in separate files
## flag=1: use name as count column name
## flag=2: use name + total count as column name
count_mut_from_tbs <- function (mut_files, name, flag=1 )
{
  all_single_mut <- NULL
  all_indel_mut <- NULL
  for(mut_file in mut_files)
  {
    tb1 <- read.table(mut_file, header=TRUE)

    single_mut <- tb1[(tb1[,"variable"] != "insertions") & (tb1[,"variable"] != "deletions") , ]

    single_mut[,"mut"] <- paste(single_mut[,"ref"],single_mut[,"variable"],sep="_")

    all_single_mut <- rbind(all_single_mut, single_mut)

    indel_mut <- tb1[(tb1[,"variable"] == "insertions") | (tb1[,"variable"] == "deletions") , ]
    all_indel_mut <- rbind(all_indel_mut, indel_mut)
  }

  ## get counts
  single_counts <- data.frame(six_mut_type_count(all_single_mut, "count"))

  ins <- sum (all_indel_mut[,"variable"] == "insertions")
  del <- sum (all_indel_mut[,"variable"] == "deletions")

  single_counts <- rbind(single_counts, c("insertions", ins))
  single_counts <- rbind(single_counts, c("deletions", del))

  single_counts[,"count"] <- as.numeric(single_counts[,"count"])

  ## count column name
  if(flag==1)
  {
    names(single_counts)[2] <- name
  }
  else if (flag==2)
  {
    names(single_counts)[2] <- paste(name,as.character(sum(single_counts[,"count"])),sep="_" )
  }

  return (single_counts)
}


## modified on 11.12.20, add indel_flag: =1: remove len(indel)>1; 0: do not remove 
## modified on 10.14.20, can specify by each sample (pass col name)
## flag=1: use strain name as the col name
## flag=2: use strain name + count as the col name
count_muts_from_tb <- function (tb, strain_name, flag=1, col_name = "strain", indel_flag=1)
{
  tb1 <- tb [tb[,col_name] == strain_name,]
  
  if(indel_flag ==1) ## only keep len(indel)==1
  {
    tb1 <- tb1[tb1[,"len"] == 1,]
  }

  single_mut <- tb1[(tb1[,"variable"] != "insertions") & (tb1[,"variable"] != "deletions") , ]

  single_mut[,"mut"] <- paste(single_mut[,"ref"],single_mut[,"variable"],sep="_")

  indel_mut <- tb1[(tb1[,"variable"] == "insertions") | (tb1[,"variable"] == "deletions") , ]

  ## count single mut spec
  single_counts <- data.frame(six_mut_type_count(single_mut, "count"))

  ## count indel
  ins <- sum (indel_mut[,"variable"] == "insertions")
  del <- sum (indel_mut[,"variable"] == "deletions")

  single_counts <- rbind(single_counts, c("insertions", ins))
  single_counts <- rbind(single_counts, c("deletions", del))

  single_counts[,"count"] <- as.numeric(single_counts[,"count"])

  ## count column name
  if(flag==1)
  {
    names(single_counts)[2] <- strain_name
  }
  else if (flag==2)
  {
    names(single_counts)[2] <- paste(strain_name, as.character(sum(single_counts[,"count"])),sep="_" )
  }

  return (single_counts)
}

## count total number of muts per sample (to compare with # of mut pooled)
count_mut_samples_from_tb0 <- function (tb)
{
  tb %>% group_by(sample) %>% summarise(sum (count)) %>% data.frame ->  sample_count

  return (sample_count)
}


count_mut_samples_from_tb <- function(mut_file)
{
  tb <- read.table(mut_file, header=TRUE)

  count_tb <- count_mut_samples_from_tb0(tb)

  return (count_tb)
}

## read from table, not from file
count_muts_from_full_tb0 <-  function (tb, strain_names, flag=1, col_name = "strain" , indel_flag=1)
{
  count_tb_all <- count_muts_from_tb(tb, strain_names[1], flag, col_name)

  if(length(strain_names) >=2 )
  {
    for(strain_name in strain_names[2:length(strain_names)])
    {
      count_tb <- count_muts_from_tb(tb, strain_name, flag, col_name,indel_flag)
      count_tb_all <- merge(count_tb_all, count_tb, by="mut")
    }
  }

  return (count_tb_all)
}

## add indel flag (remove all len(indel) > 1): 1: remove; 0: do not remove 
## count mutation spectrum from the full table
## will take in a list of all strain names
count_muts_from_full_tb <-  function (mut_file, strain_names, flag=1, col_name= "strain", indel_flag=1 )
{
  tb <- read.table(mut_file, header=TRUE)

  count_tb_all <- count_muts_from_full_tb0(tb, strain_names, flag, col_name, indel_flag)

  return (count_tb_all)
}

count_muts_from_full_tb_by_sample <- function(mut_tb, strain_name)
{
  ## extract all mutations from that strain 
  mut_tb_strain <- mut_tb [mut_tb[,"strain"] == strain_name,]
  samples <- levels(factor(mut_tb_strain[,"sample"]))
  count_tb_all <- count_muts_from_full_tb0(mut_tb_strain, samples, 1, "sample")
  return (count_tb_all)
}

## calculate pairwise euclidian square distance of mutation spectrum of different samples 
## the table is organized so that each column is a sample, each row is a mutation type 
## the first col is the mutation names
calc_mean_dis <- function(mut_tb)
{
  ncol <- ncol(mut_tb) 
  total_dis <- 0
  for(i in 2:(ncol-1) )
  {
    sample1 <- mut_tb[,i]
    for(j in (i+1): ncol)
    {
      sample2 <- mut_tb[,j]
      ## calculate pairwise euclidian distance 
      dis <- sqrt(sum((sample1-sample2)^2))
      total_dis <- total_dis + dis
    }
  }
  norm_dis <- total_dis/choose(ncol-1, 2)
  return (norm_dis)
}

## calculate mean distance for all strain 
calc_mean_dis_all_strain <- function(mut_tb)
{
  all_strains <- levels(factor(mut_tb[,"strain"]))
  dis_vec <- NULL
  for(strain in all_strains)
  {
    mut_by_sample <- count_muts_from_full_tb_by_sample(all_mut, strain)
    mut_tb_frq<-  count_2_freq_multi_col_all(mut_by_sample)
    dis <- calc_mean_dis(mut_tb_frq)
    dis_vec <- c(dis_vec, dis)
  }
  strain_dis <- data.frame(strain = all_strains, dis_vec)
  return (strain_dis)
}

