## count mutations from pooled seq data

## also will plot mut frq in histogram 

library(ggplot2)
library(gridExtra)

count_mut_tb_core <- function(mut_fluc)
{
  ## calculate mean and sd
  frq_mean <- mean(mut_fluc[,"frq"])
  frq_sd <- sd (mut_fluc[,"frq"])
  
  mut_fluc[,"count"] <- 1 ## default count = 1, unless frq > mean +2.5*sd
  
  indices <- mut_fluc[,"frq"] > frq_mean + 2.5*frq_sd
  
  multi_count_indices <- indices ## it has the final indices of mult-count mutants 
  ## indices originally contains the mult-count mutants, but later in the loop, if values have been identified, will be removed from the list (change to FALSE)
  while (! all(indices== FALSE)) 
  {
    last_indices <- indices
    
    frq_mean <- mean(mut_fluc[!multi_count_indices,"frq"])
    frq_sd <- sd (mut_fluc[!multi_count_indices,"frq"])
    
    ## recalculate indices 
    indices <- mut_fluc[,"frq"] > frq_mean + 2.5*frq_sd
    
    indices[multi_count_indices] <- FALSE
    
    multi_count_indices <- multi_count_indices | indices
  }
  
  mut_fluc[multi_count_indices, "count"] <- round (mut_fluc[multi_count_indices, "frq"] / frq_mean) 
 
  return(mut_fluc) 
}

## treat replica independently here 
## will output 2 files: 1. mut that is not present in ON file (mut_fluc) 2. mut that overlaps in ON file (mut_ON)
## 10.8 add suffix (for filter20, filter40)
count_mut_tb <- function(prefix , mut_ON_tb_file, suffix = "" , dir="./")
{
  mut_tb_file <- paste(dir, prefix, suffix, ".txt", sep="")
  mut_ON_tb_file_full <- paste(dir, mut_ON_tb_file, sep="")
  
	mut_tb <- read.table(mut_tb_file, header=TRUE)
	mut_ON_tb <- read.table(mut_ON_tb_file_full, header=TRUE)
	
	if(nrow(mut_tb)==0)
	{
	  msg <- paste("error! no mutations detected in ", mut_tb_file, sep="")
	  print (msg)
	  return (-1)
	}
	else
	{
	  mut_fluc <-  mut_tb [! (mut_tb[,"pos"] %in% mut_ON_tb[,"pos"]),]
	  
	  mut_fluc <- count_mut_tb_core(mut_fluc)
	  
	  ## then for the mut that is present in ON
	  mut_ON <- mut_tb[(mut_tb[,"pos"] %in% mut_ON_tb[,"pos"]),]
	  
	  
	  fluc_tb_name <- paste(dir, prefix, suffix, "_fluc_mut.txt", sep="")
	  on_mut_tb_name <- paste(dir, prefix, suffix, "_ON_mut.txt", sep="")
	  write.table(mut_fluc, fluc_tb_name, sep="\t",quote=FALSE, row.names=FALSE)
	  write.table(mut_ON, on_mut_tb_name, sep="\t",quote=FALSE, row.names=FALSE)
	  
	  return (0) 
	}
}



plot_pool_mut_frq <- function(file_name, title, exp_frq)
{
  
  real_mut_tb <- read.table(file_name,stringsAsFactors=FALSE, header=TRUE)
  
  
  p <- ggplot(real_mut_tb, aes(x=frq)) + geom_histogram(color="black", fill="white")
  
  p <- p + geom_vline(aes(xintercept=exp_frq),
                      color="blue", linetype="dashed", size=1) + scale_y_continuous(
                        breaks=c(1,2,3,4,5,6,7,8,9))  +labs(title=title)
  
  return (p)
  
}


count_mut_and_plot <- function(prefix , mut_ON_tb_file,  exp_frq, suffix="", dir="./")
{
    if (count_mut_tb(prefix , mut_ON_tb_file, suffix,  dir)!= -1)
    {
      fluc_tb_name <- paste(dir, prefix, suffix, "_fluc_mut.txt", sep="")
      ## plot histogram of mut 
      p <- plot_pool_mut_frq(fluc_tb_name, prefix, exp_frq)
      
      plot_name <- paste(dir, "hist/", prefix, suffix,  "_hist.pdf", sep="")
      ggsave(plot_name)
    }
}

count_mut_and_plot_all <- function(sample_info_tb, dir= "../../all_mut_wd/")
{
  ## read in the sample info  table 
  sample_info <-  read.table(sample_info_tb, stringsAsFactors = FALSE)
  
  ## only get samples which are not O.N. 
  mut_sample_tb <- sample_info [sample_info[,5] != "N",]
  
  #dir <- "../../all_mut_wd/"
  for (i in 1:nrow(mut_sample_tb))
  {
    sample_prefix <- mut_sample_tb[i,4]
    
    ## currently ON also have applied two filters, will use filter 20 file
    sample_ON <- paste(mut_sample_tb[i,5],"_filter20.txt",sep="")
    
    exp_frq <- 1/sample_info[i,2]
    
    count_mut_and_plot(sample_prefix , sample_ON, exp_frq, "_filter20", dir)
    count_mut_and_plot(sample_prefix , sample_ON, exp_frq, "_filter40", dir)
  }
}