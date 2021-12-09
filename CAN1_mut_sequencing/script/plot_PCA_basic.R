library("FactoMineR")

library(ggplot2)

library(gridExtra)

library(data.table)


library("factoextra")

library(seqinr)


rev_compl <- function(mut_tb)
{
  mut_tb[,"A_C"] <- mut_tb[,"A_C"] + mut_tb[,"T_G"]
  mut_tb[,"A_G"] <- mut_tb[,"A_G"] + mut_tb[,"T_C"]
  mut_tb[,"A_T"] <- mut_tb[,"A_T"] + mut_tb[,"T_A"]
  mut_tb[,"C_A"] <- mut_tb[,"C_A"] + mut_tb[,"G_T"]
  mut_tb[,"C_G"] <- mut_tb[,"C_G"] + mut_tb[,"G_C"]
  mut_tb[,"C_T"] <- mut_tb[,"C_T"] + mut_tb[,"G_A"]
  
  return (mut_tb[, c("A_C", "A_G","A_T", "C_A", "C_G", "C_T")])
}

## reverse complement for triplet
## each col is a mut type, each row is an individual/haplotype
rev_compl_tri <- function(mut_tb)
{
  mut_list <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
  
  ## reverse complement for single mut (use list)
  single_mut_rev_comp <- list()
  single_mut_rev_comp[["A_C"]] <- "T_G"
  single_mut_rev_comp[["A_G"]] <- "T_C"
  single_mut_rev_comp[["A_T"]] <- "T_A"
  single_mut_rev_comp[["C_A"]] <- "G_T"
  single_mut_rev_comp[["C_G"]] <- "G_C"
  single_mut_rev_comp[["C_T"]] <- "G_A"
  
  neighbor_list <- c("AA", "AC","AG","AT" , "CA","CC","CG","CT", "GA","GC", "GG","GT", "TA","TC","TG", "TT")
  
  newtb <- NULL
  col_mut <- NULL
  for (mut in mut_list)
  {
    for (neighbor in neighbor_list)
    {
      mut_triplet = paste(mut, neighbor, sep="_")
      ## reverse complement mut 
      rev_comp_nb <- c2s(rev(comp(s2c(neighbor), forceToLower = FALSE))) 
      rev_com_mut <- single_mut_rev_comp[[mut]]
      
      mut_triplet_rc <- paste(rev_com_mut,rev_comp_nb, sep="_")
      
      newcol <-  mut_tb[, mut_triplet] + mut_tb[, mut_triplet_rc] 
      col_mut <- c(col_mut, mut_triplet)
      newtb <- cbind(newtb, newcol)
    }
  }
  df <- data.frame(newtb)
  names(df) <- col_mut
  return (df)
}


## merge only 
merge_only <- function(file_vec)
{
  file <- file_vec[1]
  mut_tb_all  <- read.table(file, header=TRUE, stringsAsFactors= FALSE)
  for (file in file_vec[2:length(file_vec)])
  {
    mut_tb <- read.table(file, header=TRUE, stringsAsFactors= FALSE)
    mut_tb_all <- merge(mut_tb_all,mut_tb, by="mut" )
  }
  
  return (mut_tb_all)
}


## merge 2 or more count table to perform PCA (note that each file should have the same format)
merge_and_PCA <- function (file_vec, suffix_flag=0, rev_comp_flag=1, combine_flag=1)
{
  mut_tb_all <- merge_only(file_vec)
  
  mut.pca <- PCA_analysis_tb(mut_tb_all, suffix_flag, rev_comp_flag, combine_flag)
  return (mut.pca)
  
}




## add combine_flag =1: has 2 haplotype per individual, need to combine; =0: has one haplotype per individual, do not need to combine
## rev_comp_flag =1: single mut, needs reverse coml; =0, do not need reverse comp; =3: triplet, needs reverse compl
## this deals with a single file that each inividual has 2 counts 
## now can deal with suffix(the suffix is separated by the last "_")
PCA_analysis <- function(file, suffix_flag=0, rev_comp_flag=1, combine_flag=1)
{
  mut_tb <- read.table(file, header=TRUE, stringsAsFactors= FALSE)
  
  mut.pca <- PCA_analysis_tb(mut_tb, suffix_flag, rev_comp_flag, combine_flag)
  return (mut.pca )
}

## convert reserse complement; convert to freq 
PCA_preprocess <- function(mut_tb, suffix_flag=0, rev_comp_flag=1, combine_flag=1)
{
  ## transform the table, each line is an individual
  mut_tb_new <- t(mut_tb[,2:ncol(mut_tb)])
  
  colnames(mut_tb_new) <- mut_tb[,1]
  
  if (rev_comp_flag ==1)
  {
    mut_tb_new_rc <- rev_compl(mut_tb_new)
  } else if (rev_comp_flag ==0) ## no reverse compl
  {
    mut_tb_new_rc <- mut_tb_new
  }  else if(rev_comp_flag ==3) ## for triplet
  {
    mut_tb_new_rc <- rev_compl_tri(mut_tb_new)
  }
  
  if(combine_flag==1)
  {
    ## then add odd and even rows up 
    nows <- nrow(mut_tb_new_rc)
    odd <- mut_tb_new_rc[seq(1,nows,by= 2),]
    even <- mut_tb_new_rc[seq(2,nows,by= 2),]
    
    total <- (odd + even)/2
    rown <- rownames(odd)
  }
  else
  {
    total <- mut_tb_new_rc
  }
  
  
  if(combine_flag==1)
  {
    ## if has suffix  
    if(suffix_flag ==1)
    {
      ## find the last "_" in the strain names
      g <- regexpr("_.[^_]*$", rown)
      name_pt1 <- substr(rown, 1, g-3) ## get first part of strain name
      name_pt2 <- substr(rown, g+1, length(rown)) ## get the suffix
      strain_names <- paste(name_pt1, name_pt2, sep="+")
    }
    else
    {
      ## find the last "_" (which will have _L and _R)
      g <- regexpr("_.[^_]*$", rown)
      strain_names <-  substr(rown, 1, g-1)
    }
    rownames(total) <- strain_names
  }
  
  
  
  ## normalize the count to frequency
  rowsum <- rowSums(total)
  
  mut_tb_new_norm <- total/ rowsum
  
  return (mut_tb_new_norm)
}


## add is_scale
PCA_analysis_tb <-function(mut_tb, suffix_flag=0, rev_comp_flag=1, combine_flag=1, is_scale=FALSE)
{
  mut_tb_new_norm <- PCA_preprocess(mut_tb,suffix_flag, rev_comp_flag,combine_flag  )
  
  ## pca analysis (method 1)
  #mut.pca <- PCA(mut_tb_new_norm, scale.unit=FALSE, graph=FALSE)
  
  #PCA_tb <- data.frame(mut.pca$ind$coord)
  
  
  #ggplot(PCA_tb, aes(Dim.1,Dim.2)) + geom_point(size=1) + theme_bw() + theme( panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),  legend.title=element_blank())  
  
  #### pca using prcomp()
  mut.pca  <- prcomp( mut_tb_new_norm, scale. = is_scale)
  
  return (mut.pca)
}


###old
## include PCA analysis
#PCA_analysis <- function(file)
#{
#  mut_tb <- read.table(file, header=TRUE, stringsAsFactors= FALSE)
  
#  ## transform the table, each line is an individual
#  mut_tb_new <- t(mut_tb[,2:ncol(mut_tb)])
  
#  colnames(mut_tb_new) <- mut_tb[,1]
  
#  mut_tb_new_rc <- rev_compl(mut_tb_new)
  
  ## then add odd and even rows up 
#  nows <- nrow(mut_tb_new_rc)
#  odd <- mut_tb_new_rc[seq(1,nows,by= 2),]
#  even <- mut_tb_new_rc[seq(2,nows,by= 2),]
  
#  total <- (odd + even)/2
#  rown <- rownames(odd)
#  strain_names <- tstrsplit(rown, split="_")[[1]]
#  rownames(total) <- strain_names
  
  
  ## normalize the count to frequency
#  rowsum <- rowSums(total)
  
#  mut_tb_new_norm <- total/ rowsum
  
  ## pca analysis (method 1)
  #mut.pca <- PCA(mut_tb_new_norm, scale.unit=FALSE, graph=FALSE)
  
  #PCA_tb <- data.frame(mut.pca$ind$coord)
  
  
  #ggplot(PCA_tb, aes(Dim.1,Dim.2)) + geom_point(size=1) + theme_bw() + theme( panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),  legend.title=element_blank())  
  
  #### pca using prcomp()
  
#  mut.pca  <- prcomp(mut_tb_new_norm, scale. = FALSE)
  
#  return (mut.pca)
#}



## this function is to test two other pca functions in R
## file is the data table with individual mutation counts
PCA_vector_plot <- function(file, output_fig)
{
  
  mut.pca <- PCA_analysis (file)
  
  PCA_tb <- data.frame(mut.pca$x)
  
  ## plot PCA by individual 
  p1 <- ggplot(PCA_tb, aes(PC1,PC2)) + geom_point(size=1) + theme_bw() + theme( panel.grid.major = element_blank(),
                                                                                panel.grid.minor = element_blank(),  legend.title=element_blank())  
  
  
  #p1 <- fviz_pca_ind(mut.pca,
  #             col.ind = "#696969", # Color by the quality of representation
  #             repel = TRUE     # Avoid text overlapping
  #)
  
  ## plot vector of loadings
  p2 <- fviz_pca_var(mut.pca,
                     col.var = "contrib", # Color by contributions to the PC
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE     # Avoid text overlapping
  )
  
  ## biplot of individual and variables
  #p3 <- fviz_pca_biplot(mut.pca, repel = TRUE,
  #                col.var = "#2E9FDF", # Variables color
  #                col.ind = "#696969"  # Individuals color
  #)
  
  plots1 <-list(p1, p2)
  
  pdf(output_fig)
  
  grid.arrange(grobs = plots1, nrow=2)
  
  
  dev.off()
}




## now could contain MA spec 
## deal with the case with multiple species and mutlipe categories of mutations
## the id column in PCA file is the strainID_category
merge_PCA_tb_all <- function(pca_table, annot_tb_list, species_list, id_list)
{
  PCA_tb <- read.table(pca_table, sep=",", header=TRUE, stringsAsFactors=FALSE)
  
  ## split the id column 
  id_info <- tstrsplit(PCA_tb[,"id"],split="\\+")
  
  PCA_tb[,"newid"] <- id_info[[1]]
  PCA_tb[,"category"] <- id_info[[2]]
    
  ## read in the annotation table
  ## a lot of issues when using read.table() function
  ## https://kbroman.org/blog/2017/08/08/eof-within-quoted-string/
  
  for (i in 1:length(annot_tb_list))
  {
    annot_tb <- annot_tb_list[i]
    species <- species_list[i]
    id <- id_list [i]
    strain_annot <- read.delim(annot_tb, header=TRUE, sep=",", stringsAsFactors=FALSE)
    
    PCA_tb[PCA_tb[,"newid"] %in% strain_annot[,id],"species"] <- species
  }
  
  PCA_tb[,"pltsize"] <- 1
  
  ## for MA lines
  PCA_tb [PCA_tb[,"category"] == "MA" , "species"] <- PCA_tb [PCA_tb[,"category"] == "MA" , "newid"]
  PCA_tb [PCA_tb[,"category"] == "MA", "pltsize"] <- 4
  
  return (PCA_tb)
}


## annotate with species information
merge_PCA_tb <- function(pca_table, annot_tb_list, species_list, id_list)
{
  PCA_tb <- read.table(pca_table, sep=",", header=TRUE, stringsAsFactors=FALSE)
  
  ## read in the annotation table
  ## a lot of issues when using read.table() function
  ## https://kbroman.org/blog/2017/08/08/eof-within-quoted-string/
  
  for (i in 1:length(annot_tb_list))
  {
    annot_tb <- annot_tb_list[i]
    species <- species_list[i]
    id <- id_list [i]
    strain_annot <- read.delim(annot_tb, header=TRUE, sep=",", stringsAsFactors=FALSE)
    
    PCA_tb[PCA_tb[,1] %in% strain_annot[,id],"species"] <- species
  }
  
  
  return (PCA_tb)
}



## this function deal with one genome of individuals, but with different annotations (split the names)
annot_PCA_tb_all <-  function(pca_table, annot_tb_list, species_list, id_list)
{
  PCA_tb <- read.table(pca_table, sep=",", header=TRUE, stringsAsFactors=FALSE)
  
  ## split the id column 
  id_info <- tstrsplit(PCA_tb[,"id"],split="\\+")
  
  PCA_tb[,"newid"] <- id_info[[1]]
  PCA_tb[,"category"] <- id_info[[2]]
  
  ## read in the annotation table
  ## a lot of issues when using read.table() function
  ## https://kbroman.org/blog/2017/08/08/eof-within-quoted-string/
  
  for (i in 1:length(annot_tb_list))
  {
    annot_tb <- annot_tb_list[i]
    species <- species_list[i]
    id <- id_list [i]
    strain_annot <- read.delim(annot_tb, header=TRUE, sep=",", stringsAsFactors=FALSE)
    
    PCA_tb[PCA_tb[,"newid"] %in% strain_annot[,id],"species"] <- species
  }
  
  PCA_tb[,"pltsize"] <- 1
  

  
  return (PCA_tb)
}






## flag=0, no color, but add percent of variance explained, variable is the file name of the percent var explained 
## flag=1, color by species 
## flag=2, color by species, shape by category 
plot_PCA <- function (PCA_tb_new, title, flag, variable)
{
  
	if (flag==0)
	{
	  
	  ## read in percent of variance file 
	  percent_vars <- read.table(variable)
	  pc1_var <- percent_vars[1]
	  pc2_var <- percent_vars[2]
	  pc1_axis_str <- paste("PC1 (", round(100*pc1_var, 2), "%)", sep="")
	  pc2_axis_str <- paste("PC2 (", round(100*pc2_var, 2), "%)", sep="")
	  
	  
	## basic plot, no color
	  p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(size=1) + theme_bw() + theme( panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title) + 
	    scale_x_continuous(name=pc1_axis_str) + scale_y_continuous(name=pc2_axis_str)
  }
  else if (flag==1) 
  {
    ## discrete
    p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(size=1,aes( color=species)  ) + theme_bw() + theme( panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title)
  }
  else if (flag==2)
  {
    p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(aes( color=species, shape=category, size=pltsize)  ) + theme_bw() + theme( panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title) +
       scale_shape_manual(values = c(1,2, 13,5))
    
  }
  else if (flag==3)## also for all species, all kinds. but split category into its own plot 
  {
    p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(aes( color=species, shape=category, size=pltsize)  ) + theme_bw() + theme( panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title) +
      scale_shape_manual(values = c(1,2, 13,5)) + facet_wrap (~category, ncol = 2)  
  }
  else if (flag==4) ## same as 2, except do not scale for size 
  {
    p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(size=1, aes( color=species, shape=category)  ) + theme_bw() + theme( panel.grid.major = element_blank(),
                                                                                                                                  panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title) +
      scale_shape_manual(values = c(1,2, 13,5))
  }
  else if (flag==5) ## color by species, add percent of variance explained
  {
    ## read in percent of variance file 
    percent_vars <- read.table(variable)
    pc1_var <- percent_vars[1]
    pc2_var <- percent_vars[2]
    pc1_axis_str <- paste("PC1 (", round(100*pc1_var, 2), "%)", sep="")
    pc2_axis_str <- paste("PC2 (", round(100*pc2_var, 2), "%)", sep="")
    
    p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(size=1,aes( color=species)  ) + theme_bw() + 
      theme( panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title)  + 
      scale_x_continuous(name=pc1_axis_str) + scale_y_continuous(name=pc2_axis_str)
    
  }
  else
  {
    p <- ggplot(PCA_tb_new, aes(PC1,PC2)) + geom_point(size=1, aes( color=species, shape=category)  ) + theme_bw() + theme( panel.grid.major = element_blank(),
                                                                                                                                  panel.grid.minor = element_blank(),  legend.title=element_blank())    + labs(title = title) +
      scale_shape_manual(values = c(1,2, 13,5)) + facet_wrap (~category, ncol = 2)  
  }
  
	
	return (p)

}


## this only plots raw PCA, without color  
plot_PCA_analysis_tb_raw <- function (mut.pca,  title, pc_sign_vec)
{
  
  PCA_tb <- data.frame(mut.pca$x)
  
  eigs <- mut.pca$sdev^2
  percent_vars = eigs / sum(eigs)
  
  PCA_tb[,"ID"] <- rownames(PCA_tb)
  
  
  pc1_axis_str <- paste("PC1 (", round(100*percent_vars[1], 2), "%)", sep="")
  pc2_axis_str <- paste("PC2 (", round(100*percent_vars[2], 2), "%)", sep="")
  
  
  p <- ggplot(PCA_tb, aes(pc_sign_vec[1]*PC1,pc_sign_vec[2]*PC2)) + 
    geom_point(size=1) + theme_bw() + theme( panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank(),  legend.title=element_blank()) +
    labs(title = title) + scale_x_continuous(name=pc1_axis_str) + scale_y_continuous(name=pc2_axis_str)
  
  
  return (p)
  
  
}

## add option to reverse x or y axis 
## reverse_vec = c(-1,1) means reverse x axis but not reverse y axis  
plot_PCA_vec <- function(mut.pca, reverse_vec)
{
  p <- fviz_pca_var(mut.pca,
                    col.var = "contrib", # Color by contributions to the PC
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE     # Avoid text overlapping
  )
  p <- p  + theme( panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
  
  if (reverse_vec[1] == -1)
  {
    p <- p + scale_x_reverse()
  }
  
  if (reverse_vec[2] == -1)
  {
    p <- p + scale_y_reverse()
  }
  return (p)
}

