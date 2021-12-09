## plot the basic Scer mut PCA

library(dplyr)
library(ggrepel)
library(cowplot)
library(gridExtra)

source("plot_PCA_basic.R")
source("choose_strains_basic.R")

## this part only run once 
## first need to subsample strains from
## each pop, no more than 30
## and then perform PCA 
## and project to the rest 


plot_PCA_analysis <- function( file, title, rev_comp_flag, annot_file, combine_flag ,pc_sign_vec, color_flag=1)
{
  mut.pca <- PCA_analysis (file, 0, rev_comp_flag, combine_flag )
  p <- plot_PCA_analysis_tb(mut.pca, title, rev_comp_flag, annot_file, combine_flag ,pc_sign_vec, color_flag)
  return (p)
}



## color_flag=1: only color African beer, TW, French dairy, wine pop; =2: color a few more pop
## file is the mut count
## plot PCA analysis (done in R)
## pc_sign_vec is a vector that indicates the sign to plot PC1 and PC2. e.g. if vec is -1, -1, will plot -PC1, -PC2.
plot_PCA_analysis_tb <- function (mut.pca, PCA_tb,  title, annot_file ,pc_sign_vec, color_flag=1)
{
  
  eigs <- mut.pca$sdev^2
  percent_vars = eigs / sum(eigs)
  
  PCA_tb[,"ID"] <- rownames(PCA_tb)
  
  ## read in annot tb 
  annot <- read.table(annot_file, sep=",", header=TRUE)
  
  ## merge annot 
  PCA_tb  <- merge(PCA_tb, annot, by="ID")
  
  PCA_tb[,"new_group"] <- "Other"
  
  ## add an extra column for the presence and abscence African beer 
  if(color_flag==1)
  {
    
    PCA_tb [PCA_tb[,"group"] == 1 , "new_group" ] <- "1.Wine/European" 
    PCA_tb [PCA_tb[,"group"] == 5 , "new_group" ] <- "5.French dairy" # 
    PCA_tb [PCA_tb[,"group"] == 6 , "new_group" ] <- "6.African beer" # African_beer
    PCA_tb [PCA_tb[,"group"] == 17 , "new_group" ] <- "17.Taiwanese" # Taiwanese
  }
  else if(color_flag==2)
  {
    PCA_tb [PCA_tb[,"group"] == 1 , "new_group" ] <- "1.Wine/European" 
    PCA_tb [PCA_tb[,"group"] == 5 , "new_group" ] <- "5.French dairy" # 
    PCA_tb [PCA_tb[,"group"] == 6 , "new_group" ] <- "6.African beer" # African_beer
    PCA_tb [PCA_tb[,"group"] == 17 , "new_group" ] <- "17.Taiwanese" # Taiwanese
    PCA_tb [PCA_tb[,"group"] == 3 , "new_group" ] <- "3.Brazilian Bioethanol" 
    PCA_tb [PCA_tb[,"group"] == 7 , "new_group" ] <- "7.Mosaic beer" # African_beer
  

    PCA_tb [PCA_tb[,"group"] == 13 , "new_group" ] <- "13.African palm wine" #
    
  
    PCA_tb [PCA_tb[,"group"] == 25 , "new_group" ] <- "25. Sake"
    
    ## add a few more groups 
    PCA_tb [PCA_tb[,"group"] == 24 , "new_group" ] <- "24. Asian Islands"
    
    PCA_tb [PCA_tb[,"group"] == 26, "new_group" ] <- "26. Asian Fermentation"
    
    PCA_tb [PCA_tb[,"group"] == 8, "new_group" ] <- "8. Mixed origin"
    
  }
  else if(color_flag==3)
  {
    PCA_tb [PCA_tb[,"group"] == 13 , "new_group" ] <- "13.African palm wine" #
    PCA_tb [PCA_tb[,"group"] == 15 , "new_group" ] <- "15.CHN"
   
  }
  
  
  color_vec <- c( "1.Wine/European"= "#0BC346", 
                  "3.Brazilian Bioethanol"="#CC6677",
                  "5.French dairy"= "#88CCEE", 
                  "6.African beer"=  "#000000", 
                  "7.Mosaic beer"="#2F30EE",
                  "8. Mixed origin"= "#FFEE08",
                  "13.African palm wine"= "#D55E00",
                  "17.Taiwanese" =  "#DC0101",
                  "24. Asian Islands"= "#AA4499",
                  "25. Sake" = "#332288",
                  "26. Asian Fermentation"= "#CAB23C",
                  "Other"= "#b3b3b3"
  )
  
  PCA_tb[,"new_group"] <- factor(PCA_tb[,"new_group"], 
                                 levels= c("1.Wine/European", "3.Brazilian Bioethanol",
                                            "5.French dairy", "6.African beer", "7.Mosaic beer",
                                           "8. Mixed origin", "13.African palm wine", "17.Taiwanese",
                                           "24. Asian Islands", "25. Sake",
                                            "26. Asian Fermentation", 
                                            "Other" ))  
  
  ## remove numbers 
  new_names <- c("Wine/European", "Brazilian Bioethanol", 
                 "French dairy", "African beer", "Mosaic beer", "Mixed origin", 
                 "African palm wine", "Taiwanese", "Asian Islands", 
                 "Sake",
                 "Asian Fermentation", 
                 "Other" ) 
  
  
  ## get mean and sd for PC1 and PC2 
  ## will plot strain IDs if they are mean+-1.8sd away
  mean1 <- mean(PCA_tb[,"PC1"])
  sd1 <- sd(PCA_tb[,"PC1"])
  mean2 <- mean(PCA_tb[,"PC2"])
  sd2 <- sd(PCA_tb[,"PC2"])

  idx1 <- (PCA_tb[,"PC1"] <= mean1 - 1.8*sd1 ) | (PCA_tb[,"PC1"] >= mean1 + 1.8*sd1 )
  idx2 <- (PCA_tb[,"PC2"] <= mean2 - 1.8*sd2 ) | (PCA_tb[,"PC2"] >= mean2 + 1.8*sd2 )
  
  ## generate a label column only to label the outliers 
  PCA_tb[,"label"] <- ""
  
  PCA_tb[idx1 | idx2,"label"] <- PCA_tb[idx1 | idx2, "ID"]
  
  pc1_axis_str <- paste("PC1 (", round(100*percent_vars[1], 2), "%)", sep="")
  pc2_axis_str <- paste("PC2 (", round(100*percent_vars[2], 2), "%)", sep="")
  
  
  p <- ggplot(PCA_tb, aes(pc_sign_vec[1]*PC1, pc_sign_vec[2]*PC2,  label=label)) + 
    geom_hline(yintercept=0, linetype="dashed", color="grey" ) + ## plot dotted lines first, then plot dots
    geom_vline(xintercept=0, linetype="dashed", color="grey" ) + 
    geom_point(size=1.8, colour="#5a5a5a", shape=21) + 
    aes(fill= new_group) +
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), 
                        legend.title=element_blank(), axis.title=element_text(size=12,face="bold"),
                        axis.text = element_text(size=11, face="bold")) +
    scale_fill_manual(values= color_vec, labels =new_names )  + 
    labs(title = title) + scale_x_continuous(name=pc1_axis_str) + scale_y_continuous(name=pc2_axis_str, limits=c(-0.19, 0.35) ) +
    geom_text_repel(size=3, colour = "black") 
  
  return (p)
}




## merge the strain count table with the annotation table
## and only use strains in the annotation table
annot_tb <- read.table ("../data/strain_subset_final_AC4_anno_no_close.csv", header =TRUE, sep=",",
                        stringsAsFactors = FALSE)

all_count <- read.table("../data/Scer_single_no_close_AC4_indi_mut_count_strain.txt", header =TRUE)
all_count_mut <- all_count[-c(13,14), ]

## exclude strains which have total mutation count <8
count_tb <- all_count_mut[,2:ncol(all_count_mut)]
total <- colSums(count_tb)
all_count_no_outlier <- count_tb[,total >=8 ]
outlier <- names(count_tb[,total <8 ])
all_count_no_outlier_tb <- cbind("mut"=all_count_mut[,1] ,all_count_no_outlier )


## exclude outlier from annot tb 
subset <- annot_tb [! (annot_tb[,"ID"] %in% outlier), "ID"]
mut_count_tb <- all_count_no_outlier_tb[,c("mut", subset)]

## to be consistent with common var 
## exclude introgressed strains 
introg <- annot_tb[annot_tb[,"group"] ==2 | annot_tb[,"group"] ==9,"ID"]

subset2 <- subset [!subset %in% introg]
mut_count_tb2 <- mut_count_tb[,c("mut", subset2)]

 
## will output strains used in this PCA 
## exclude introgressed strains 
## previous version: Scer_AC4_PCA_strains.txt (every thing else the same, except including the introgressed strains)
write.table(names(mut_count_tb2), "../data/Scer_AC4_PCA_strains_v2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


setEPS()
postscript("../plot/fig1_rare_var_PCA_raw.eps", width=5.9, height=6.2)


mut.pca <- PCA_analysis_tb(mut_count_tb2,0, 1,0, FALSE)
  
  
  
## project PCA on a smaller dataset to the whole strains 
mut_tb_new_norm <- PCA_preprocess(mut_count_tb2,0, 1,0)
PCA_tb_new <- scale(mut_tb_new_norm, mut.pca$center, mut.pca$scale) %*% mut.pca$rotation %>% data.frame
  
  
single_p <- plot_PCA_analysis_tb(mut.pca, PCA_tb_new,"", 
                                   "/net/harris/vol1/project/yeast_1002/data/strain_anno_group.txt",c(-1,-1),2)
  
  
## plot vector 
p2 <- plot_PCA_vec(mut.pca,c(-1,-1))
  
plots <-list( single_p,  p2)
  

  
#https://stackoverflow.com/questions/36198451/specify-widths-and-heights-of-plots-with-grid-arrange
plot_grid(single_p, p2,  align = "v", nrow = 2, rel_heights = c(2/3, 1/3))


dev.off()

