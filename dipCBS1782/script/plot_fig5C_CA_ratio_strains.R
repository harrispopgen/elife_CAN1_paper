source("/net/harris/vol1/home/pyjiang/mut_R/mut_basic_funcs.R")
library(ggplot2)
library(ggrepel)
library(cowplot)

calc_mut_from_tb <- function(tb)
{
  
  tb[,"count"] <- 1
  
  ## calculate mut counts
  mut_count <- data.frame(six_mut_type_count(tb,"count"))
  
  return (mut_count)
}


calc_mut_from_file <- function(file, colnames )
{
	tb <-read.table(file,sep=",", stringsAsFactors=FALSE)

	names(tb) <- colnames
	mut_count<- calc_mut_from_tb(tb)
	
	return (mut_count)
}


## get freq from table 
calc_mut_freq_from_file <- function(file, colnames )
{
  tb <-read.table(file,sep=",", stringsAsFactors=FALSE)
  
  names(tb) <- colnames
  tb[,"count"] <- 1
  
  ## calculate mut freq
  mut_freq <- data.frame(six_type_mut_freq_only(tb,"count"))
  
  return (mut_freq)
}

## get count or frequency of a specific mutation type and the other type 
## return type: 2 element vector 
## flag=0: count; flag=1: freq
## input: a mutation table count with 6 types of mutations  
count_mut_one_type_frq <- function(mut_tb, mut_type, flag)
{
  ## count mutations of specific type 
  mut_type_count <- mut_tb[mut_tb[,"mut"] == mut_type, "count"] 
  ## count mutations of all 
  mut_all_count <- sum(mut_tb[, "count"] )
  
  rest <- mut_all_count-mut_type_count
  
  if(flag==0)
  {
    return (c(mut_type_count, rest))
  }
  else
  {
    return (c(mut_type_count/mut_all_count, rest/mut_all_count))
  }
}


## get raw counts 
get_all_count_tb <- function(prefix, list_of_strains, list_of_AC )
{
  all_mut_tb <- NULL
  colnames = c("chr", "pos", "mut", "AF", "AC", "strains")
  
  for(i in 1:length(list_of_strains))
  {
    for (j in 1:length(list_of_AC))
    {
      strain <- list_of_strains[i]
      AC <- list_of_AC[j]
      file <- paste(prefix, strain, "_DAC_", toString(AC), "_single_mut.csv" ,sep="")
      mut_tb <- calc_mut_from_file(file, colnames)
      
      mut_tb[,"strain"] <- strain
      mut_tb[,"AC"] <- AC
      
      all_mut_tb <- rbind(all_mut_tb, mut_tb)
    }
  }
  return (all_mut_tb)
}


## get frq of one type of mutation 
## modify prefix to the mutation AC prefix (for all mut)
get_frq_one_type_mut <- function(prefix, list_of_strains, list_of_AC, mut_type)
{
  all_mut_frq_tb <- NULL
  colnames = c("chr", "pos", "mut", "AF", "AC", "strains")
  
  for (j in 1:length(list_of_AC))
  {
    AC <- list_of_AC[j]
    mut_AC_file <- paste(prefix, toString(AC),".csv", sep="")
    mut_tb_full <- read.table(mut_AC_file, sep=",", header=TRUE) 
    names(mut_tb_full) <- c("chr","pos","mut","AC","sample")
    for(i in 1:length(list_of_strains))
    {
      strain <- list_of_strains[i]
      #file <- paste(prefix, strain, "_DAC_", toString(AC), "_single_mut.csv" ,sep="")
      subtb <- mut_tb_full [mut_tb_full[,"sample"] == strain,]
      mut_tb <- calc_mut_from_tb(subtb)
      mut_frq <- count_mut_one_type_frq(mut_tb, mut_type, 1)[1]
      total_mut <- sum(mut_tb[,"count"])
      
      vec <- data.frame(ratio = mut_frq,strain = strain, AC= AC, total = total_mut)
      
      all_mut_frq_tb <- rbind(all_mut_frq_tb, vec)
    }
  }
 
  return (all_mut_frq_tb)
  
}


## get freq 
get_all_freq_tb <- function(prefix, list_of_strains, list_of_AC )
{
  all_mut_tb <- NULL
  colnames = c("chr", "pos", "mut", "AF", "AC", "strains")
  
  for(i in 1:length(list_of_strains))
  {
    for (j in 1:length(list_of_AC))
    {
      strain <- list_of_strains[i]
      AC <- list_of_AC[j]
      file <- paste(prefix, strain, "_DAC_", toString(AC), "_single_mut.csv" ,sep="")
      mut_tb <- calc_mut_freq_from_file(file, colnames)
      
      mut_tb[,"strain"] <- strain
      mut_tb[,"AC"] <- AC
      
      all_mut_tb <- rbind(all_mut_tb, mut_tb)
    }
  }
  return (all_mut_tb)
}


## plot mutation count/frq comparision 
plot_mut_bar <- function(prefix, list_of_strains, list_of_AC ,figname )
{
  color_vec <- c("C_A"= "#E69F00", "A_C"= "#999999", "A_T" ="#999999", "A_G"="#999999",
                 "C_T"= "#999999", "C_G" = "#999999")
  frq_tb <- get_all_freq_tb(prefix, list_of_strains, list_of_AC  )
  
  mut_tb <- get_all_count_tb(prefix, list_of_strains, list_of_AC)
  mut_count_AC <- mut_tb %>% group_by(strain, AC) %>% summarise(sum(count)) %>% data.frame
  names(mut_count_AC) <- c("strain", "AC", "count")
  ## assign counts to C_G mutation (so it only shows up on top of that bar)
  mut_count_AC[,"mut"] <- "C_G"
  mut_count_AC[,"display"] <- paste("n=", mut_count_AC[,"count"],sep="")
  
  frq_tb_w_count <- merge(frq_tb, mut_count_AC, by=c("strain", "AC", "mut"), all=TRUE)
  ## will 
  
  p <- ggplot(frq_tb_w_count, aes(x=mut, y=freq, fill=mut)) + geom_bar(stat="identity") + 
    facet_grid(vars(strain), vars(AC) ) +
    scale_fill_manual(values = color_vec) + 
    geom_text(aes(label=display) , vjust=-8, hjust=0.6, size=3)
  ggsave(figname)
}

## for strains with potential elevation in C->A mut 
#plot_mut_bar("../data/beer_mut_", c("AEA", "YJM271", "AEQ", "AAR", "SP007", "SP008"), c(4,8,12,16) ,"../plot/beer_mut_comp1.png")

## for strains that do not show elevation in C->A mut 

#plot_mut_bar("../data/beer_mut_", c("AAQ", "BelleSaison","YJM193" , "AQH"), c(4,8,12,16) ,"../plot/beer_mut_comp2.png")


## negative: AAQ, BelleSaison, YJM193



## plot ratio of C->A mutations with different levels of AC (allele count in the total population)
## add total count for each strain as well 
plot_mut_CA_ratio_count <- function(prefix, list_of_strains, list_of_AC ,figname, color_vec, wid, heigh )
{
  C_A_frq_tb <- get_frq_one_type_mut(prefix, list_of_strains, list_of_AC, "C_A")
  
  ## note: need to specify "color" to make legend (if specify group variable, will not show legend!)
  ## then will manually change color 
  p <- ggplot(C_A_frq_tb, aes(x=AC, y=ratio, color = strain, label=total)) + 
    geom_line(lty="longdash", size=1.1) + 
    geom_point() +
    labs(y="C>A ratio", x= "Allele Count") + scale_color_manual(values= color_vec)
  
  p2 <- p + geom_text_repel() +
    theme_bw(base_size = 11) + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),  legend.title=element_blank(),
                        axis.title=element_text(size=12,face="bold"),
                        axis.text = element_text(size=11, face="bold") ) 
  ggsave(figname, width=wid, height=heigh)
}

## do not plot the total counts of mutations on the dot plot, but 
## just output to a table 
plot_mut_CA_ratio_count_v2 <- function(prefix, list_of_strains, list_of_AC ,figname, color_vec,  heigh, ratio )
{
  C_A_frq_tb <- get_frq_one_type_mut(prefix, list_of_strains, list_of_AC, "C_A")
  
  ## note: need to specify "color" to make legend (if specify group variable, will not show legend!)
  ## then will manually change color 
  p <- ggplot(C_A_frq_tb, aes(x=AC, y=ratio, color = strain)) + 
    geom_line(lty="longdash", size=1.1) + 
    geom_point() +
    labs(y="C>A ratio", x= "Allele Count") + scale_color_manual(values= color_vec)
  
  p2 <- p + 
    theme_bw(base_size = 11) + theme( panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),  legend.title=element_blank(),
                                      axis.title=element_text(size=12,face="bold"),
                                      axis.text = element_text(size=11, face="bold") ) 
  
  #ggsave(figname, width=wid, height=heigh)
  
  ## convert to short tb 
  #C_A_frq_tb_short <- dcast(C_A_frq_tb, strain ~ AC)
  
  strain_show <- c("AAR", "AEQ", "SACE_YAG", "BRM", "CBS1782")
  C_A_frq_tb_sub <- C_A_frq_tb[C_A_frq_tb[,"strain"] %in% strain_show, ]
  
  p3 <- ggplot(C_A_frq_tb_sub, aes(x=AC, y=strain, color = strain, label= total)) +
    geom_tile(aes(fill= strain), colour = "grey50" )+
    geom_text(col = "black") +
    theme_bw() +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position = "none") + scale_fill_manual(values= color_vec)
  
  #ggsave("../plot/fig5_mut_comp_CA_ratio_count.eps", width=wid, height=1.7)
  
  p_out <- plot_grid(p3, p2,  align = "v", nrow = 2, rel_heights = c(0.25, 0.75),  axis = 'lr')
  #save_plot("../plot/fig1_bar.svg", p_out, base_height=5, base_asp= 0.7)
  
  save_plot(figname, p_out, base_height=heigh, base_asp= ratio)
}



## for y1002 strains 

color_vec_y1002 <- c("BRM"= "#446455"  , "SACE_YAG"=  "#E69F00", "AEQ"=  "#0072B2", "AAR" =  "#C93312" , "CBS1782" = "#7f2d80","AAQ"= "#999999",
                      "SACE_YAB" = "#999999" , "AQH"= "#E69F77" , "AFP"="#999999" , "AFA"= "#999999", "AFB"="#999999")

plot_mut_CA_ratio_count_v2 ("../data/mut_info_AC_", c("BRM", "SACE_YAG", "AEQ", "AAR", "SACE_YAB", "CBS1782", "AQH", "AFP", "AFA", "AFB", "AAQ" ), 
                         c(2,4,6,8) ,"../plot/fig5B_mut_comp_CA_ratio_add_CBS1782.svg", color_vec_y1002, 5, 1.05)






