## plot fig1 of the Scer paper
## include suppl fig 1) PCA for syn mutations 
## 2) mutation spectrum differences from all natural variants versus MA, and singletons. [previous 1B]


########
## mutation PCA
#######


library(data.table)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggrepel)


source("/net/harris/vol1/home/pyjiang/mut_R/plot_PCA_basic.R")
source("/net/harris/vol1/home/pyjiang/mut_R/mut_basic_funcs.R")
source("/net/harris/vol1/home/pyjiang/mut_R/call_expr_mut_basic.R")


## 1. will calculate PCA based on a subset of strains (<=30); if project with all, too ugly (mostly are wine strains)
## 2. add loadings 

### not use 
PCA_project_to_all <- function(subset_file, title, rev_comp_flag, annot_file, combine_flag, pc_sign_vec, color_flag, all_indi_file)
{
  mut.pca <- PCA_analysis (subset_file, 0, rev_comp_flag, combine_flag )
  
  mut_count_tb <- read.table(all_indi_file, header=TRUE, stringsAsFactors = FALSE)
  mut_tb_new_norm <- PCA_preprocess(mut_count_tb,0, rev_comp_flag, combine_flag)
  
  PCA_tb_new <- scale(mut_tb_new_norm, mut.pca$center, mut.pca$scale) %*% mut.pca$rotation %>% data.frame
  p <- plot_PCA_analysis_tb(mut.pca, PCA_tb_new, title, rev_comp_flag, annot_file, combine_flag ,pc_sign_vec, color_flag)
  return (p)
  
}


plot_PCA_analysis <- function( file, title, rev_comp_flag, annot_file, combine_flag ,pc_sign_vec, color_flag=1)
{
  mut.pca <- PCA_analysis (file, 0, rev_comp_flag, combine_flag )
  PCA_tb <- PCA_tb <- data.frame(mut.pca$x)
  p <- plot_PCA_analysis_tb(mut.pca, PCA_tb, title, rev_comp_flag, annot_file, combine_flag ,pc_sign_vec, color_flag)
  return (p)
}

## modified to input with a new PCA_tb 
## color_flag=1: only color African beer, TW, French dairy, wine pop; =2: color a few more pop
## file is the mut count
## plot PCA analysis (done in R)
## pc_sign_vec is a vector that indicates the sign to plot PC1 and PC2. e.g. if vec is -1, -1, will plot -PC1, -PC2.
plot_PCA_analysis_tb <- function (mut.pca,  PCA_tb, title, rev_comp_flag, annot_file, combine_flag ,pc_sign_vec, color_flag=1, plot_outlier_flag=FALSE)
{
  
  eigs <- mut.pca$sdev^2
  percent_vars = eigs / sum(eigs)
  
  PCA_tb[,"ID"] <- rownames(PCA_tb)
  
  ## read in annot tb 
  annot <- read.table(annot_file, sep=",", header=TRUE)
  
  #PCA_tb  <- cbind(PCA_tb, annot )
  
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
    
    
    #PCA_tb [PCA_tb[,"group"] == 12 , "new_group" ] <- "12.West African cocoa" #
    
    PCA_tb [PCA_tb[,"group"] == 13 , "new_group" ] <- "13.African palm wine" #
    
    
   # PCA_tb [PCA_tb[,"group"] == 14 , "new_group" ] <- "14-18.CHN" ## note: not 15!
  #  PCA_tb [PCA_tb[,"group"] == 16 , "new_group" ] <- "14-18.CHN"
    
   # PCA_tb [PCA_tb[,"group"] == 15 , "new_group" ] <- "14-18.CHN"
  #  PCA_tb [PCA_tb[,"group"] == 18 , "new_group" ] <- "14-18.CHN"
    
    
    #PCA_tb [PCA_tb[,"group"] == 23 , "new_group" ] <- "23. North American Oak"
    
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
                  "17.Taiwanese" =  "#DC0101",
                  "5.French dairy"= "#88CCEE", 
                  "6.African beer"=  "#000000", 
                  "7.Mosaic beer"="#2F30EE",
                  "8. Mixed origin"= "#FFEE08",
                  "3.Brazilian Bioethanol"="#CC6677",
                  #"14-18.CHN" = "#b3b3b3",
                  "13.African palm wine"= "#D55E00",
                  "26. Asian Fermentation"= "#CAB23C",
                  "25. Sake" = "#332288",
                  "24. Asian Islands"= "#AA4499",
                  #"12.West African cocoa" = "#FAD510",
                  #"23. North American Oak"="#3CAEA3",
                  "Other"= "#b3b3b3"
                  )
  
  PCA_tb[,"new_group"] <- factor(PCA_tb[,"new_group"], levels= c("1.Wine/European", "17.Taiwanese",
			"5.French dairy", "6.African beer", "7.Mosaic beer",
			"8. Mixed origin",
			"3.Brazilian Bioethanol",
			#"14-18.CHN", 
			"13.African palm wine", "26. Asian Fermentation", "25. Sake",
			"24. Asian Islands", 
			#"12.West African cocoa", "23. North American Oak",
			"Other" ))  
  
  ## remove numbers 
  new_names <- c("Wine/European", "Taiwanese",
                 "French dairy", "African beer", "Mosaic beer", "Mixed origin", 
                 "Brazilian Bioethanol",
                 #"14-18.CHN", 
                 "African palm wine", "Asian Fermentation", "Sake",
                 "Asian Islands", 
                 #"12.West African cocoa", "23. North American Oak",
                 "Other" ) 
  
  pc1_axis_str <- paste("PC1 (", round(100*percent_vars[1], 2), "%)", sep="")
  pc2_axis_str <- paste("PC2 (", round(100*percent_vars[2], 2), "%)", sep="")
  
  ## to be consistent to fig6 
  
  if(plot_outlier_flag==FALSE)
  {
    p <- ggplot(PCA_tb, aes(pc_sign_vec[1]*PC1,pc_sign_vec[2]*PC2)) + 
      geom_hline(yintercept=0, linetype="dashed", color="grey" ) + ## plot dotted lines first, then plot dots
      geom_vline(xintercept=0, linetype="dashed", color="grey" ) + 
      geom_point(size=1.8, colour="#5a5a5a", shape=21) + 
      aes(fill= new_group) + 
      theme_bw(base_size = 11) + theme( panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),  
                                        legend.title=element_blank(),
                                        axis.title=element_text(size=12,face="bold"),
                                        axis.text = element_text(size=11, face="bold") )+
      scale_fill_manual(values= color_vec, labels =new_names )  + 
      labs(title = title) + scale_x_continuous(name=pc1_axis_str) + scale_y_continuous(name=pc2_axis_str)
  }
  else
  {
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
      labs(title = title) + scale_x_continuous(name=pc1_axis_str) + scale_y_continuous(name=pc2_axis_str, limits=c(-0.2, 0.35) ) +
      geom_text_repel(size=3, colour = "black") 
  }
  

  
  return (p)
  #ggsave(out_name, width=6,height=4)
  
}

setEPS()
postscript("../plot/fig1_PCA_new.eps", width=5.3, height=5.2)

#svg("../plot/fig1_PCA.svg", width=5.3, height=5.2)

#mut.pca <- PCA_analysis ("../../data/yeast_single_indi_mut_count.txt", 0, 1, 1 )

## read in the count for all individuals 
## then get the subset of strains that overlap with rare variants PCA (excluding introgressed lines)

## read in count for all strains 
## then only get the subset that is used for AC4. 
all_mut_count <- read.table("../data/Scer_single_no_close_indi_mut_count_strain.txt", header=TRUE)


## read in the strains used for AC4 
strains <- read.table("../data/strain_subset_final_AC4_no_close_raw.txt")

strains_L <- paste(strains[,1],"L", sep="_")
strains_R <- paste(strains[,1],"R", sep="_")
strains_full <- c(rbind(strains_L, strains_R))

## get individuals only in the strains list specified  
mut_count <- all_mut_count[,c("mut",strains_full )]

mut.pca <- PCA_analysis_tb (mut_count, 0, 1, 1 )
 PCA_tb <- data.frame(mut.pca$x)

single_p <- plot_PCA_analysis_tb(mut.pca,PCA_tb, "",1, 
                              "../../data/strain_anno_group.txt",1,c(1,-1),2)

p2 <- plot_PCA_vec(mut.pca, c(1,-1))

plot_grid(single_p, p2,  align = "v", nrow = 2, rel_heights = c(2/3, 1/3))

dev.off()

#ggsave("../plot/fig1_PCA.svg", width=7, height=5)


###################
## PCA plot from syn variants (for suppl)
###################

pdf("../plot/figS1_PCA_syn.pdf", width=5.3, height=5.2)

mut_syn_all <- read.table("../data/Scer_single_no_close_indi_syn_mut_count.txt", header=TRUE)
syn_mut_count <- mut_syn_all[,c("mut",strains_full )]
  
mut.pca <- PCA_analysis_tb (syn_mut_count, 0, 1, 1 )
PCA_tb <- data.frame(mut.pca$x)

single_p <- plot_PCA_analysis_tb(mut.pca,PCA_tb, "",1, 
                                 "../../data/strain_anno_group.txt",1,c(-1,-1),2)

p2 <- plot_PCA_vec(mut.pca, c(-1,-1))

title <- ggdraw() + draw_label("Synonymous variants")


plot_grid(title, single_p, p2,  align = "v", ncol=1, rel_heights = c(0.08, 0.62, 0.3))

dev.off()



###################
## PCA plot from randomly drawn variants from decendents (for suppl)
###################

pdf("../plot/figS2_PCA_rand.pdf", width=5.3, height=5.2)

mut_rand_all <- read.table("../data/Scer_single_no_close_indi_mut_count_rand.txt", header=TRUE)
rand_mut_count <- mut_rand_all[,c("mut",strains_full )]

mut.pca <- PCA_analysis_tb (rand_mut_count, 0, 1, 1 )

PCA_tb <- data.frame(mut.pca$x)

single_p <- plot_PCA_analysis_tb(mut.pca,PCA_tb, "",1, 
                                 "../../data/strain_anno_group.txt",1,c(1,1),2)

p2 <- plot_PCA_vec(mut.pca, c(1,1))

title <- ggdraw() + draw_label("Randomly assign mutation to a decendent")

plot_grid(title, single_p, p2,  align = "v", ncol=1, rel_heights = c(0.08, 0.62, 0.3))

dev.off()


##################
## PCA plot for singletons (suppl)
##################
## singleton tb for haploids
singleton_tb1 <- read.table("/net/harris/vol1/project/yeast_1002/extremely_rare/homozygous/haploid/singleton_tb.txt", header =TRUE)
## singleton tb for homozygous diploids
singleton_tb2 <- read.table("/net/harris/vol1/project/yeast_1002/extremely_rare/homozygous/diploid_homo/singleton_tb.txt", header =TRUE)
## singleton tb for heterozygous (regardless of ploidy)
singleton_tb3 <- read.table("/net/harris/vol1/project/yeast_1002/extremely_rare/heterozygous/singleton_tb.txt", header =TRUE)

## combine singleton tb 
singleton_tb <- rbind(singleton_tb1, singleton_tb2, singleton_tb3)

### read in the strains used in rare variants (AC=4) PCA
#strains <- read.table("../data/Scer_AC4_PCA_strains.txt", header=TRUE)

## only include strains that are used in AC=4 PCA for comparison 
singleton_tb_sub <- singleton_tb [singleton_tb[,"ind"] %in% strains[,1],]



singleton_tb_sub[,"mut"] <- paste(singleton_tb_sub[,"ref"], singleton_tb_sub[,"alt"], sep="_")
names(singleton_tb_sub)[5] <- "strain"
singleton_tb_sub[,"count"] <- 1

count_per_strain <- six_mut_type_count_strains(singleton_tb_sub, "count")

## also exclude strains which have total mutation count <8
count_tb <- count_per_strain[,2:ncol(count_per_strain)]
total <- colSums(count_tb)
all_count_no_outlier <- count_tb[,total >=8 ]
outlier <- names(count_tb[,total <8 ])
all_count_no_outlier_tb <- cbind("mut"=count_per_strain[,1] ,all_count_no_outlier )

mut.pca2 <-  PCA_analysis_tb (all_count_no_outlier_tb, 0,0,0)
PCA_tb2 <- data.frame(mut.pca2$x)

pdf("../plot/figS3_PCA_singleton.pdf", width=6, height=5.9)

single_p <- plot_PCA_analysis_tb(mut.pca2,  PCA_tb2, "",1, 
                                 "../../data/strain_anno_group.txt",1,c(1,1),2, TRUE)

p2 <- plot_PCA_vec(mut.pca2, c(1,1))

title <- ggdraw() + draw_label("Singleton")

plot_grid(title, single_p, p2,  align = "v", ncol=1, rel_heights = c(0.08, 0.67, 0.25))

dev.off()


##################
## PCA plot for triplet (suppl)
##################


mut_triplet_all <- read.table("../data/Scer_triplet_no_close_indi_mut_count_strain.txt", header=TRUE)
triplet_mut_count <- mut_triplet_all[,c("mut",strains_full )]

mut.pca <- PCA_analysis_tb (triplet_mut_count, 0, 3, 1 )

PCA_tb <- data.frame(mut.pca$x)

triplet_p <- plot_PCA_analysis_tb(mut.pca,PCA_tb, "",1, 
                                 "../../data/strain_anno_group.txt",1,c(-1,1),2) + ggtitle("Triplet context")


ggsave("../plot/figS_PCA_triplet.pdf", width=5.3, height=4)



###################
## bar plot comparing MA and natural spectrum
## with singletons (fig S5)
###################

## read in natural var 
Scer_rand_mut_count <- read.table("../../data/yeast_single_indi_mut_count_all_rand.txt", header =TRUE)

Scer_nat_counts <- count_all_mut_cols (Scer_rand_mut_count, "All natural") ## counts 

Scer_nat <- count_all_mut_2_frq(Scer_rand_mut_count,"All natural") ## frq


## read in counts for young singletons and all singletons from haploids/homo diploids or other strains

Scer_homo_singleton_count <- read.table("../../data/homo_singleton_mut_count.txt", header =TRUE)

Scer_other_singleton_count <- read.table("../../data/other_singleton_mut_count.txt", header =TRUE)


Scer_homo_singleton_young <- count_2_freq_multi_col(Scer_homo_singleton_count, "young_singleton","Young singletons")
Scer_homo_singleton_all <- count_2_freq_multi_col(Scer_homo_singleton_count, "all_singleton","All singletons")

Scer_other_singleton_young <- count_2_freq_multi_col(Scer_other_singleton_count, "young_singleton","young singleton other")
Scer_other_singleton_all <- count_2_freq_multi_col(Scer_other_singleton_count, "all_singleton","all singleton other")



## mutation freq of CAN1 gene from previous literature (summarized in Lynch 2008)
mut_list <- c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T")
can1_frq <- c(0.066, 0.094, 0.075, 0.269, 0.137, 0.358)

can1_frq_tb <- data.frame(mut = mut_list, CAN1 = can1_frq)


## record mutation freq from Zhu et al 2014 (fig 3)
Zhu_frq <- c(0.11, 0.144, 0.063, 0.182, 0.152, 0.35)
Zhu_frq_tb <- data.frame(mut = mut_list, CAN1 = Zhu_frq)


## read in Sharp mut data
de_novo_mut_tb <- read.table("/net/harris/vol1/project/yeast_1002/data/Sharp_MA_replic.txt",header=TRUE)

## use haploid spectrum (use all hap)

hap_mut_tb <- de_novo_mut_tb [ (de_novo_mut_tb[,"Ploidy"] == "hap") & (de_novo_mut_tb[,"RDH54"] == "wt") ,]


hap_mut_tb[,"mut"] <- paste(hap_mut_tb[,"Ref"], hap_mut_tb[,"Alt"], sep="_")

hap_mut_tb[,"count"] <- 1

hap_mut_frq_tb <- six_type_mut_freq_only(hap_mut_tb,"count")

names(hap_mut_frq_tb) <- c("mut", "MA" )

## get counts
hap_mut_count_tb <- six_mut_type_count(hap_mut_tb,"count")
names(hap_mut_count_tb) <- c("mut", "MA")



## combine frq tables
full_mut_tb <- merge(Scer_nat,Scer_homo_singleton_all, by="mut" )
full_mut_tb <- merge(full_mut_tb,Scer_homo_singleton_young, by="mut" )

#full_mut_tb <- merge(full_mut_tb,Scer_other_singleton_all, by="mut" )
#full_mut_tb <- merge(full_mut_tb,Scer_other_singleton_young, by="mut" )

full_mut_tb <- merge(full_mut_tb, hap_mut_frq_tb, by= "mut")
#full_mut_tb <- merge(full_mut_tb, can1_frq_tb , by="mut")



## generate two plots, one include late replic MA, the other does not include 

full_mut_tb_long <- reshape2::melt(full_mut_tb)

## calculate transition to transversion ratio 

ts_tb <- full_mut_tb_long %>% filter(mut=="A_G" | mut=="C_T" )
tv_tb <- full_mut_tb_long %>% filter(mut=="A_C" | mut=="A_T" |  mut=="C_A" |  mut=="C_G" )

ts_sum <- ts_tb %>% group_by(variable)   %>% summarise(sum(value))
names(ts_sum) <- c("variable", "ts")
tv_sum <- tv_tb %>% group_by(variable)   %>% summarise(sum(value))
names(tv_sum) <- c("variable", "tv")

ts_tv <- merge(ts_sum, tv_sum, by= "variable")
ts_tv[,"ts_tv"] <- ts_tv[,"ts"]/ ts_tv[,"tv"]

color_vec <- c("A_C" ="#1674bb" , "A_G" = "#f2a3f0", "A_T" = "#37343f",
               "C_A"="#4daf4a", "C_G" ="#85e1de", "C_T" = "#c4250c")

color_vec2 <- c("A_C" ="#1674bb" , "A_G" = "#f2a3f0", "A_T" = "#37343f",
               "C_A"="#4daf4a", "C_G" ="#85e1de", "C_T" = "#c4250c", "ts_tv" = "#F98400")

ts_tv_sub <- ts_tv[,c("variable", "ts_tv")]

#ts_tv_sub_long <- melt(ts_tv_sub)
#names(ts_tv_sub_long) <- c("variable", "mut", "value")
#ts_tv_sub_long2 <- ts_tv_sub_long[,c("mut", "variable","value" )]

#full_mut_tb_long2 <- rbind(full_mut_tb_long, ts_tv_sub_long2)

full_mut_tb_long2 <- merge(full_mut_tb_long,ts_tv_sub, by= "variable" )


## plot ratio and ts/tv with lines
## two y axis 
## https://www.r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
p1 <- ggplot(data=full_mut_tb_long2, aes(x=variable)) +
      geom_point(aes( y=value, colour=mut)) +
      geom_line( aes( y=value, colour=mut, group=mut), linetype="dashed") +
      geom_point(aes(y=ts_tv/10), color="#F98400", shape=1) + 
      geom_line (aes(y=ts_tv/10,  group=1), color="#F98400", linetype=1) +
      scale_color_manual( values= color_vec2) + 
      scale_y_continuous(
        # Features of the first axis
        name = "Frequency",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~. *10, name="ts/tv")
      ) + theme_bw() + 
     theme( panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
            axis.title.y.right = element_text(color = "#F98400"),
            axis.title=element_text(size=12,face="bold"),
            axis.text = element_text(size=12, face="bold")) +
      labs(x="")


## barplot
p2 <- ggplot(data=full_mut_tb_long, aes(x=variable, y=value, fill=mut))+ geom_bar(position= "fill", stat="identity") +
  scale_fill_manual( values= color_vec) +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 90),
         axis.title=element_text(size=12,face="bold"),
         axis.text = element_text(size=12, face="bold")) +
  labs(y= "Frequency", x="", fill="Mutation")


## use cowplot to align two plots together (not including legends)
## https://stackoverflow.com/questions/16255579/how-can-i-make-consistent-width-plots-in-ggplot-with-legends/16258375#16258375
## https://www.rdocumentation.org/packages/cowplot/versions/1.0.0/topics/save_plot
p_out <- plot_grid(p1, p2,  align = "v", nrow = 2, rel_heights = c(0.3, 0.7))

save_plot("../plot/figS5_bar.eps", p_out, base_height=5, base_asp= 0.7)



##############
# plot mutation spectra comparing CAN1 (from Lang) and MA
# 6.16.21, add MA data from Zhu et al
##############

## read in mut counts from my CAN1 data
can1_all_mut <- read.table("/net/harris/vol1/project/yeast_1002/CAN1_mut_sequencing/data/mut_count_tb_no_indel.txt", header=TRUE)

##  count for GIL and Lang
can1_all_frq<-  count_2_freq_multi_col_all(can1_all_mut)
GIL_Lang_can1 <- can1_all_frq[,c("mut", "GIL_311", "Lang")]

can1_vs_MA1 <- merge(GIL_Lang_can1, hap_mut_frq_tb, by="mut") ## Sharp
can1_vs_MA2 <- merge(can1_vs_MA1, Zhu_frq_tb, by="mut") ##Zhu

names(can1_vs_MA2) <- c("mut", "CAN1 (this study)", "CAN1 (Lang et al.)","MA (Sharp et al.)", "MA (Zhu et al.)")

can1_vs_MA_long <- reshape2::melt(can1_vs_MA2)

## bar plot

p <- ggplot(data=can1_vs_MA_long, aes(x=mut, y=value, fill=mut))+ geom_bar(stat="identity") +
  scale_fill_manual( values= color_vec) +facet_wrap(~variable, ncol=2 ) +
  scale_x_discrete(labels=c("insertions"= "ins", "deletions" = "del")) +
  theme_bw() +
  theme( strip.placement = "outside",
         strip.background = element_blank(),
         strip.text = element_text( size = 12 ),
         axis.line = element_line(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 90),
         axis.title=element_text(size=12,face="bold"),
         axis.text = element_text(size=12, face="bold")) + 
  labs(y= "Frequency", x="", fill="Mutation")

ggsave("../plot/figS_CAN1_MA.pdf", width=5.5, height=3.5)



