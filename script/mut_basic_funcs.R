## basic functions to calculate mutations

library(dplyr)
library(reshape2)
library(seqinr)

## change the reverse complement to mut which begin with A and C
cat_rev_comp <- function (full_tb)
{
	## then change the reverse complement mut to be the same (only keep those begin with A and C)
	full_tb [full_tb["mut"] =="T_G","mut"] <- "A_C"
	full_tb [full_tb["mut"] =="T_C","mut"] <- "A_G"
	full_tb [full_tb["mut"] =="T_A","mut"] <- "A_T"
	full_tb [full_tb["mut"] =="G_T","mut"] <- "C_A"
	full_tb [full_tb["mut"] =="G_C","mut"] <- "C_G"
	full_tb [full_tb["mut"] =="G_A","mut"] <- "C_T"

	return (full_tb)
}



## revere complement function for triplet 
## require to have 3 columns 1) ref 2) mut 3) five_prime base 4) three_prime base
## note: ref means the base before mutation, not necessarily ref base
cat_rev_comp_triplet <- function (full_tb)
{
  full_tb[,"neighbor"] <- paste(full_tb[,"five_prime"], full_tb[,"three_prime"],sep="")
  
  ## change mut column
  full_tb <- cat_rev_comp(full_tb)
  
  full_tb[,"ori_neighbor"] <- full_tb[,"neighbor"]
  
  bases = c("A","C","G","T")
  
  ## for those that are converted rev compl 
  for (mid_base in c("G","T"))
  {
    idx <- full_tb["ref"] == mid_base
    
    for (x in bases)
    {
      for (y in bases)
      {
        neighbor <- paste(x,y,sep="")
        
        #print (neighbor)
        rev_comp_nb <- c2s(rev(comp(s2c(neighbor), forceToLower = FALSE))) 
        
        context <- paste(x,mid_base, y, sep="")
        rev_comp_ct <- c2s(rev(comp(s2c(context), forceToLower = FALSE)) )
        
        idx2 <- full_tb["ori_neighbor"] == neighbor
        full_tb[idx & idx2,"neighbor"] = rev_comp_nb
        full_tb[idx & idx2,"context"] = rev_comp_ct
        
      }
    }
  }
  ## for those which are the original (will add context)
  for (mid_base in c("A","C"))
  {
    idx <- full_tb["ref"] == mid_base
    
    full_tb[idx,"context"] <- paste(full_tb[idx,"five_prime"], full_tb[idx,"ref"],full_tb[idx,"three_prime"],sep="" )
  }
  
  full_tb["triplet"] <- paste(full_tb[,"mut"], full_tb[,"neighbor"],sep="_")
  
  return (full_tb)
  
  
}


## create an empty vector with 96 mut type (triplet)
empty_triplet <- function()
{
  mut_list <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
  neighbor_list <- c("AA", "AC","AG","AT" , "CA","CC","CG","CT", "GA","GC", "GG","GT", "TA","TC","TG", "TT")
  
  triplet_vec <- NULL 
  for (mut in mut_list)
  {
    for (neighbor in neighbor_list)
    {
      mut_triplet = paste(mut, neighbor, sep="_")
      triplet_vec <- c(triplet_vec,mut_triplet )
    }
  }
  return (triplet_vec)
}

## empty triplet with another column of pop name 
empty_triplet_pop <- function(pop_name)
{
  triplet <- empty_triplet()
  df <- data.frame(pop=pop_name, triplet = triplet)
  return (df)
}



## convert mutation types  (3 letters(context)_mut) to 2 digit(do not consider nearby site)
## mut_col is the column number or name of the 3 digit mutation
## as a result, add a column named "mut" into the full_tb
six_mut_type <- function(full_tb, mut_col)
{

	mut_type_vec <- full_tb[, mut_col]

	before <- substr(as.character(mut_type_vec),2,2)
	after <-  substr(as.character(mut_type_vec),5,5)
	#full_tb ["before"] <- before
	#full_tb ["after"] <- after

	## concatenate before and after together
	full_tb["mut"] <- paste(before,after,sep="_")
	
	full_tb  <- cat_rev_comp(full_tb)
	
	return (full_tb)
}


## the full_tb has a column with "mut", A_C means A mutates to C
six_mut_type_count <- function (full_tb, count_col_name)
{
	
	full_tb  <- cat_rev_comp(full_tb)

	mut_tb <- full_tb [,c(count_col_name,"mut")]

	## calculate the number of mutations in each mutation category

	mut_tb %>% group_by(mut) %>% summarise_all(sum) -> mut_count
	
	## if there is zero count for a mutation type, will add 0 to it 
	if(nrow(mut_count) <6)
	{
		mut_list <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
		for (mut in mut_list)
		{
			if(nrow( mut_count[ mut_count[,"mut"]==mut, ])==0 )
			{
				
				mut_count  %>% add_row( mut = mut, count = 0) -> mut_count
			}
		}
	}

	return (mut_count)
}

## new function, will calculate mut count by strains 
## input has a column with all strain names (column name: strain)
## will return a matrix, with strain names in columns, mutation type in rows.
six_mut_type_count_strains <- function (full_tb, count_col_name)
{
  full_tb  <- cat_rev_comp(full_tb)
  
  mut_tb <- full_tb [,c(count_col_name,"mut","strain")]
  
  
  ## calculate the number of mutations in each mutation category
  
  mut_tb %>% group_by(mut, strain) %>% summarise_all(sum) %>%  data.frame -> mut_count_tb 
  
  ## replace NA with 0
  mut_count_matrix <-  reshape2::dcast(mut_count_tb, mut~strain)
  mut_count_matrix[is.na(mut_count_matrix)] <- 0
  
  return (mut_count_matrix)
  
}



#################
## note: this function is wrong!
#################
## count the number of mutations for all 96 types 
## must have field: mut (single mutation type, A_C for example), neighbor (5' and 3' bases: for example AT: 5' A, 3' T)
mut_type_triplet_count <- function(full_tb, count_col_name)
{
	full_tb  <- cat_rev_comp(full_tb)
		
	full_tb["triplet"] <- paste(full_tb[,"mut"], full_tb[,"neighbor"],sep="_")
	
	mut_tb <- full_tb [,c(count_col_name, "triplet")]

	
	mut_tb %>% group_by(triplet) %>% summarise_all(sum) -> mut_triplet_count
	
	## if there is zero count for a mutation type, will add 0 to it 
	if(nrow(mut_triplet_count) < 96)
	{
		mut_list <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
		neighbor_list <- c("AA", "AC","AG","AT" , "CA","CC","CG","CT", "GA","GC", "GG","GT", "TA","TC","TG", "TT")
		
		for (mut in mut_list)
		{
			for (neighbor in neighbor_list)
			{	
				mut_triplet = paste(mut, neighbor, sep="_")
				if(nrow( mut_triplet_count[ mut_triplet_count[,"triplet"]== mut_triplet, ])==0 )
				{
					mut_triplet_count  %>% add_row( triplet = mut_triplet, count = 0) -> mut_triplet_count
				}
			}
		}
	}
	
	return (mut_triplet_count)
}
 
## six mut count per chr  
 


## calculate freq of six mut type (the table will have both count and freq)
six_type_mut_freq <- function (full_tb, count_col_name)
{
	mut_count <- six_mut_type_count(full_tb, count_col_name)
	mut_count[,"freq"] <- mut_count[, count_col_name]/sum(mut_count[, count_col_name])
	return (mut_count)
}


## calculate freq of six mut type but only return with freq (drop count)
six_type_mut_freq_only <- function (full_tb, count_col_name)
{
	mut_count <- six_mut_type_count(full_tb, count_col_name)
	mut_count[,"freq"] <- mut_count[, count_col_name]/sum(mut_count[, count_col_name])
	mut_count[,count_col_name] <- NULL
	return (mut_count)
}

## convert count to freq 
## if flag=1, use count_col_name as output. if flag=0, do not .
six_type_mut_count2frq <- function(mut_count, count_col_name, flag=1)
{
  mut_count[,"freq"] <- mut_count[, count_col_name]/sum(mut_count[, count_col_name])
  mut_count[,count_col_name] <- NULL
  
  if (flag==1)
  {
    names(mut_count)[2] <- count_col_name
  }
  
  return (mut_count)
}



## mut_type_vec is a vector with mutation types (3 letters(context)_mut)
## count_col_name is the column name which has the count of each type of 3 letter mutation types
six_mut_type_count_3base <-function(mut_col, full_tb, count_col_name)
{
	
	full_tb <- six_mut_type(full_tb, mut_col)
	
	mut_count <- six_mut_type_count(full_tb, count_col_name)
	
	return (mut_count)
}


## calculate transition to transversion ratio 
## full_tb has a column of "mut", which has the six type mutation 
## will add a column of ts_tv, indicating 
ts_tv <- function(full_tb, group_col)
{
	full_tb[,"ts"] <-0
	full_tb[,"tv"] <-0
	 
	full_tb[full_tb["mut"] == "A_G","ts"] <-1
	full_tb[full_tb["mut"] == "C_T","ts"] <-1
	full_tb[full_tb["mut"] == "A_C","tv"] <-1
	full_tb[full_tb["mut"] == "A_T","tv"] <-1
	full_tb[full_tb["mut"] == "C_G","tv"] <-1
	full_tb[full_tb["mut"] == "C_A","tv"] <-1

	## note: if want to pass argument to dplyr functions, add "_" to the function! (such as group_by here)
	full_tb  %>% group_by_(group_col) %>%  summarise(sum(ts),sum(tv)) -> ts_tv_tb
	return (ts_tv_tb)
}


## function to calculate the total mutation count, as well as for each mutation type, given a data frame with annotation
count_mut_from_anno <- function(mut_tb)
{

	## convert 3 digit mutation to one (add column mut)
	#mut_tb <- six_mut_type(mut_tb,3)

	mut_tb["count"] <-1 ## add count to make it easier to sum the desired type

	## convert to shot format
	mut_tb_mut_type <- dcast(mut_tb, chr+ pos + gene_ID +count ~ mut, value.var="count")


	## replace NA with 0
	mut_tb_mut_type[is.na(mut_tb_mut_type)] <- 0


	## then get count by gene name
	## sum total mutations per gene; sum for each of six type of mut

    mut_tb_mut_type %>% group_by(gene_ID) %>% summarise(sum(count),sum(A_C),sum(A_G),sum(A_T),sum(C_A),sum(C_G),sum(C_T))   -> count_by_gene

	names(count_by_gene) <- c("gene_ID","mut_count","A_C","A_G","A_T","C_A","C_G","C_T")
	
	return (count_by_gene)
}


## this deals with two columns (one is mutation type, and the other is the count)
count_2_freq <- function(mut_count, count_col_name, new_name)
{
  mut_count[,"freq"] <- mut_count[, count_col_name]/sum(mut_count[, count_col_name])
  mut_count[,count_col_name] <- NULL
  names(mut_count)[2] <- new_name
  return (mut_count)
}


## 12.16.19 add flag to deal with triplet as well, context_flag=1: single; context_flag=3, triplet. 
## 12.13.19 (change Mut_type to mut)
## from now on, will use mut in all cases
## just count all mut by combining all column counts
## add flag=0: count all columns; flag=1: exclude insertions and deletions
count_all_mut_cols <- function(mut_tb, name, flag=0, context_flag=1)
{
  
  if(flag==1)
  {
    mut_tb1 <- subset(mut_tb,  mut!= "insertions" ) 
    mut_tb <- subset(mut_tb1,  mut!= "deletions" ) 
  }
  
  mut_count_tb <- rowSums(mut_tb[,2:ncol(mut_tb)])
  mut_type <- mut_tb[,1]
  mut_count_tb_new <- data.frame(mut=mut_type, count = mut_count_tb, stringsAsFactors = FALSE )
  if(context_flag ==1)
  {
    mut_count_tb_new2 <- as.data.frame (six_mut_type_count(mut_count_tb_new, "count"))
  }
  else if (context_flag ==3)
  {
    mut_triplet <- mut_count_tb_new[,"mut"]
    mut_count_tb_new[,"mut"] <- substr(mut_triplet,1,3)
    mut_count_tb_new[,"neighbor"] <- substr(mut_triplet,5,6)
    mut_count_tb_new2 <- as.data.frame (mut_type_triplet_count(mut_count_tb_new, "count"))
    names(mut_count_tb_new2)[1] <- "mut"
  }
  
  names(mut_count_tb_new2)[2] <- name
  mut_count_tb_new2[,1] <- as.character(mut_count_tb_new2[,1])
  
  return (mut_count_tb_new2)
}


## this mut_tb file contains mutations for each individual. it will output the counts of all mutations for each mutation type combined
count_all_mut_2_frq <- function(mut_tb, name, flag=0)
{
  mut_count_tb_new2 <- count_all_mut_cols(mut_tb, name,flag)
  
  mut_frq_tb <- count_2_freq(mut_count_tb_new2, name, name)
  return (mut_frq_tb)
}





## this mut_tb file contains multiple columns of different mutation counts (such as different species, or different categories, such as singletons, etc)
## will convert mutation counts to frequency by specify specific
count_2_freq_multi_col <- function(mut_count, count_col_name, new_name)
{
  subtb <- mut_count[,c("mut",count_col_name)]
  mut_freq <- count_2_freq(subtb,count_col_name,new_name )
  return (mut_freq)
}

## this mut_tb file contains multiple columns of different mutation counts
## first column is "mut", then from second column to the end are different sample names.
count_2_freq_multi_col_all <- function (mut_count)
{
  samples <- ncol (mut_count)  
  
  col_name <- names(mut_count)[2]
  subtb <- mut_count[,c("mut",col_name)]
  count_2_freq(subtb,col_name,col_name )
  mut_freq_all <- count_2_freq(subtb,col_name,col_name )
  
  for (i in 3:samples)
  { 
      col_name <- names(mut_count)[i]
      subtb <- mut_count[,c("mut",col_name)]
      count_2_freq(subtb,col_name,col_name )
      mut_freq <- count_2_freq(subtb,col_name,col_name )
      mut_freq_all <- merge(mut_freq_all, mut_freq, by="mut")
  }
  
  return (mut_freq_all)
}

## mut is the row names, rather than the first column 
count_2_freq_multi_col_all2 <- function (mut_count)
{
  col_sum <- colSums(mut_count)
  
  mut_frq <- sweep(mut_count,2,col_sum,"/")

  return (mut_frq)
}

