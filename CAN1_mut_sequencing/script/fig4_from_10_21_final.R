library(ggplot2)
library(cowplot)

library(gridExtra)
library(grid)
library(ggrepel)



source("call_expr_mut_basic.R")
source("plot_PCA_basic.R")

## update on 12.30.20, use the mutation file that has been corrected for position locations 


#############
# print out total number of mutations per category: MNM, single mut, indels 
#############

## calculate the occurrence of mutations for each strain
## record sites that have multiple mutations (either from the same pool or different pool)


## note that when couting mut spec for each strain, only include 1bp indels. 
tb <- read.table("../data/all_mut_combined_coord_converted.txt", header=TRUE)

## count the #of >1bp indels, 1bp indels, single nucleotide mutations 

single_mut <- tb[(tb[,"variable"] != "insertions") & (tb[,"variable"] != "deletions") , ]
indel_mut <- tb[(tb[,"variable"] == "insertions") | (tb[,"variable"] == "deletions") , ]
indel_mut_1bp <- indel_mut [indel_mut[,"len"] == 1,]

## count total #of single mut (4561)
single_sub_mut_count <- sum(single_mut[,"count"])

## count total #of indels (1bp)
indel_mut_1bp_count <- sum(indel_mut_1bp[,"count"])
indel_mut_count <- sum(indel_mut[,"count"])
indel_mut_count_long <- indel_mut_count-indel_mut_1bp_count

## read in MNM files (for the total count)
double_MNM <- read.table("../MNM/final_MNM_double.txt", header=TRUE)
complex_MNM <- read.table("../MNM/final_MNM_complex.txt", header =TRUE)

MNM_count <- nrow(double_MNM) + nrow(complex_MNM)

all_mut_count <- single_sub_mut_count + indel_mut_count + MNM_count

# print 
print(paste("single substitution count:", single_sub_mut_count))
print(paste("total indel count:", indel_mut_count))
print(paste("indel 1bp count:", indel_mut_1bp_count))
print(paste("indel >1bp count:", indel_mut_count_long))
print(paste("MNM count:", MNM_count))
print(paste("all mut count:", all_mut_count))


## compare single mutants called in my experiments versus those in Lang et al. (only consider single bp substitutions) 
## how many are now sites and how many are old sites 
Lang_mut <- read.table("../data/Lang_CAN1_mut.txt", header =TRUE)
Lang_SNV <- Lang_mut[(Lang_mut[,"from"] != "-") &  (Lang_mut[,"to"] != "-"),]
#Lang_SNV[,"mut"] <- paste(Lang_SNV[,"pos"], Lang_SNV[,"from"], Lang_SNV[,"to"], sep="_")
#single_mut[,"mut"] <- paste(single_mut[,"pos"], single_mut[,"ref"], single_mut[,"variable"],sep="_")

Lang_SNV_count <- sum(Lang_SNV[,"count"])
print(paste("total number of SNV in Lang:", Lang_SNV_count))

## exclude one mutation which occurs at the insertion seq, not in ref. 
single_mut2 <- single_mut[single_mut[,"special"] !="yes",]

unique_sites <- unique(single_mut2[,"new_pos"])
unique_sites_prev <- unique(Lang_SNV[,"pos"])

print (paste("total number of unique sites in my dataset:", length( unique_sites)))
print (paste("total number of unique sites in Lang's dataset:", length(unique_sites_prev)))

## check overlap sites 
print (paste("number of overalp sites:", sum(unique_sites_prev %in% unique_sites)))

## new sites discovered in my study 
print (paste("number of unique sites in my study:", sum(! (unique_sites %in% unique_sites_prev))))

## find positions that are uniquely in Lang's study 
only_in_Lang <- unique_sites_prev[! (unique_sites_prev %in% unique_sites),]

###############
# perform chi-square tests on mutation spectrum, also plot 4B. and supplementary figs
###############


mut_count_tb <- count_muts_from_full_tb0(tb,
                                        c("alpha54", "AGM", "AAA", "ACS", "AHH", "AHC", "AEF","ADN","ADF",
                                          "AEQ", "AAR", "ACM", "CEV", "AGR", "GIL","AAB","AAI", "AFH"), 2)


lang_mut_tb <- read.table("../data/CAN1_Lang_paper_mut_count.txt", header=TRUE)

mut_count_tb <- merge(mut_count_tb, lang_mut_tb)

## write it to files (to send to Alan)
## write the full count table 
write.table(mut_count_tb, "../data/mut_count_tb.txt",quote=FALSE,sep="\t", row.names=FALSE)
write.table(mut_count_tb[1:6,], "../data/mut_count_tb_no_indel.txt",quote=FALSE,sep="\t", row.names=FALSE)


## chi-sq test against alpha54
## p=0.05 cutoff: 
## not significant: AAA, ACS, AHC, ACM, AGR, GIL
## significant but mild: AGM, AHH, AEF, ADN, ADF, CEV, AAB, Lang
## highly significant (p<10-12): AEQ, AAR

## and record p values into mut_count_tb
cols <- names(mut_count_tb)

## create a new table that records p value as well as significant label "*"
sig_tb<- data.frame(matrix(nrow=4, ncol=length(cols)))
colnames(sig_tb) <- cols
rownames(sig_tb) <- c("pval_w_indel_label", "pval_w_indel", "pval_no_indel_label", "pval_no_indel")

for (i in 1:18)
{
  print (cols[i+2])
  tb <- data.frame(x = mut_count_tb[,"alpha54_316"], mut_count_tb[,i+2])
  chisq_tb <- chisq.test(tb)
  print (chisq_tb)
  sig_tb["pval_w_indel",i+2] <- chisq_tb$p.value
  if(chisq_tb$p.value >= 0.05)
  {
    sig_tb["pval_w_indel_label",i+2] <- "NS"
  }
  else if(chisq_tb$p.value >= 0.001)
  {
    sig_tb["pval_w_indel_label",i+2] <- "*"
  }
  else if(chisq_tb$p.value >= 0.0001)
  {
    sig_tb["pval_w_indel_label",i+2] <- "**"
  }
  else
  {
    sig_tb["pval_w_indel_label",i+2] <- "***"
  }
}


## if only consider SNPs, no indels
## not significant: AGM, AAA, ACS, AHC, ADF, ACM, AGR, GIL, Lang
## significant but mild: AHH, AEF, ADN, CEV, AAB
## highly significant (p<10-10): AEQ, AAR

mut_count_tb_snp <- mut_count_tb[1:6,]
for (i in 1:18)
{
  print (cols[i+2])
  tb <- data.frame(x = mut_count_tb_snp[,"alpha54_316"], mut_count_tb_snp[,i+2])
  chisq_tb <- chisq.test(tb)
  print (chisq_tb)
  
  sig_tb["pval_no_indel",i+2] <- chisq_tb$p.value
  if(chisq_tb$p.value >= 0.05)
  {
    sig_tb["pval_no_indel_label",i+2] <- "NS"
  }
  else if(chisq_tb$p.value >= 0.001)
  {
    sig_tb["pval_no_indel_label",i+2] <- "*"
  }
  else if(chisq_tb$p.value >= 0.0001)
  {
    sig_tb["pval_no_indel_label",i+2] <- "**"
  }
  else
  {
    sig_tb["pval_no_indel_label",i+2] <- "***"
  }
}

## only get the label 
sig_tb_keep <-  sig_tb[c("pval_w_indel_label", "pval_no_indel_label"), ]
sig_tb_keep[,"mut"] <- NULL
sig_tb_keepv <- data.frame(t(sig_tb_keep))
sig_tb_keepv[,"variable"] <- rownames(sig_tb_keepv)

sig_tb_keepv_indel <- sig_tb_keepv[,c("variable", "pval_w_indel_label")]
sig_tb_keepv_no_indel <- sig_tb_keepv[,c("variable","pval_no_indel_label")]




## plot frequency

mut_tb_frq<-  count_2_freq_multi_col_all(mut_count_tb)

## plot freq
mut_tb_sub_frq_long <- reshape2::melt(mut_tb_frq)




## calculate enrichment of C>A mutations in AAR and AEQ
# fold changes of C>A mutations to C>T mutations in AAR
0.64890282/0.19435737
## in AEQ
0.605577689/0.167330677

## if using C>T as the baseline, to estimate the overall elevation of mutation due to excessive C>A
## in AAR, = 0.7
(0.64890282/0.19435737 -1) * 0.29746835
## in AEQ = 0.78
(0.605577689/0.167330677 -1) * 0.29746835


############
# suppl fig S8-S9.
############

######### will keep the current order, and just place different stars to represent 
######### p value significance for these 
## change the orders of plots
## first show the controls, then non-significan strains, followed by significant ones
strain_order <- c("alpha54_316", "GIL_311", "Lang", # CONTROL
                  "AAA_294", "ACS_339", "AHC_303", "ACM_240", "AGR_253", "AAI_327","AAB_231",
                  "AGM_330", "AHH_346",  "AEF_328","ADN_320",  "ADF_292", "CEV_236",  "AFH_254",
                  "AEQ_252", "AAR_316") 

color_vec <- c("A_C" ="#1674bb" , "A_G" = "#f2a3f0", "A_T" = "#37343f",
               "C_A"="#4daf4a", "C_G" ="#85e1de", "C_T" = "#c4250c",
               "insertions" ="#E8E09C" , "deletions" = "#EAA92B")

mut_tb_sub_frq_long[,"variable"] <- factor(mut_tb_sub_frq_long[,"variable"], levels= strain_order)

mut_tb_sub_frq_long_indel_sig <- merge(mut_tb_sub_frq_long, sig_tb_keepv_indel, by="variable")

p <- ggplot(data=mut_tb_sub_frq_long_indel_sig, aes(x=mut, y=value, fill=mut, label=pval_w_indel_label))+ geom_bar(stat="identity") +
  geom_text(vjust = 0) +
  scale_fill_manual( values= color_vec) +facet_wrap(~variable, ncol=4 ) +
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
         axis.text = element_text(size=12, face="bold"))

#ggsave("../plot/mut_freq_new_10_21_20_new.pdf", width=12, height=8)

ggsave("../plot/figS_mut_spec_bar_w_indel.pdf", width=12, height=8)

## snp only 
mut_tb_frq_snp<-  count_2_freq_multi_col_all(mut_count_tb_snp)

mut_tb_sub_frq_snp_long <- reshape2::melt(mut_tb_frq_snp)

mut_tb_sub_frq_snp_long[,"variable"] <- factor(mut_tb_sub_frq_snp_long[,"variable"], levels= strain_order)

mut_tb_sub_frq_long_no_indel_sig <- merge(mut_tb_sub_frq_snp_long, sig_tb_keepv_no_indel, by="variable")


p <- ggplot(data=mut_tb_sub_frq_long_no_indel_sig, aes(x=mut, y=value, fill=mut, label=pval_no_indel_label))+ geom_bar(stat="identity") +
  geom_text(vjust = 0) +
  scale_fill_manual( values= color_vec) +facet_wrap(~variable, ncol=4 ) +
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
         axis.text = element_text(size=12, face="bold"))

#ggsave("../plot/mut_freq_new_10_21_20_new.pdf", width=12, height=8)

ggsave("../plot/figS_mut_spec_bar_no_indel.pdf", width=12, height=8)


###########
# suppl fig S10
###########

######### generate mutation spectrum of two batches of AAR (showing relatively high and low rates) from the experiment
all_mut_raw <- read.table("../all_mut_wd_excl_MNM/all_mut_combined.txt", header=TRUE)
## select only AAR
AAR_mut <- all_mut_raw[all_mut_raw[,"strain"] == "AAR",]

## separate AAR muts (of higher or lower rate estimates)
## higher: from r5, AAR_r23_2 (8.15.20) lower: from others 
r5 <- AAR_mut[,"sample"] == "AAR_r5_1A" | AAR_mut[,"sample"] == "AAR_r5_2"
high2 <- AAR_mut[,"sample"] == "AAR_r23_2"
high <- r5 | high2
AAR_high <- AAR_mut [ high,]
AAR_low <- AAR_mut [ !high,]

AAR_high_count <- count_muts_from_tb(AAR_high, "AAR", 2)
AAR_low_count <- count_muts_from_tb(AAR_low, "AAR", 2)

AAR_tb <- merge(AAR_high_count,AAR_low_count, by="mut" )

## high and low rates, respectively
AAR_tb2 <- AAR_tb[,c("AAR_83", "AAR_233")]
names(AAR_tb) <- c("mut", "high rate (n=83)", "low rate (n=233)")

print (chisq.test(AAR_tb2))

print (chisq.test(AAR_tb2[1:6,])) 

AAR_mut_tb_frq<-  count_2_freq_multi_col_all(AAR_tb)

AAR_mut_tb_frq_long <- reshape2::melt(AAR_mut_tb_frq)

#color_vec <- c("A_C" ="#ff7f00" , "A_G" = "#377eb8", "A_T" = "#984ea3", "C_A"="#4daf4a", "C_G" ="#ffff33", "C_T" = "#e41a1c",
#               "insertions" ="#ED8D8B" , "deletions" = "#9D4F29")


color_vec <- c("A_C" ="#1674bb" , "A_G" = "#f2a3f0", "A_T" = "#37343f",
               "C_A"="#4daf4a", "C_G" ="#85e1de", "C_T" = "#c4250c",
               "insertions" ="#E8E09C" , "deletions" = "#EAA92B")

p <- ggplot(data=AAR_mut_tb_frq_long, aes(x=mut, y=value, fill=mut))+ geom_bar(stat="identity") +
  scale_fill_manual( values= color_vec) +facet_wrap(~variable, ncol=4 ) +
  scale_x_discrete(labels=c("insertions"= "ins", "deletions" = "del")) + 
  ggtitle("AAR mutation spectrum with high or low rates")


ggsave("../plot/fig_S8_AAR_high_low.pdf", width=6.5, height=3.5)



###############
## fig 4C
###############
## perform a PCA on the mutation frq of strains 

mut.pca <- PCA_analysis_tb(mut_count_tb, rev_comp_flag=0 ,combine_flag=0)

## plot PCA result 
PCA_tb <- data.frame(mut.pca$x)

eigs <- mut.pca$sdev^2
percent_vars = eigs / sum(eigs)

PCA_tb[,"ID"] <- tstrsplit(rownames(PCA_tb), split="_")[[1]]

## read in annot tb 
annot <- read.table("/net/harris/vol1/project/mut_survey/species/Scer/data/strain_anno_group.txt", sep=",", header=TRUE)

#PCA_tb  <- cbind(PCA_tb, annot )

## merge annot 
PCA_tb  <- merge(PCA_tb, annot, by="ID", all.x=TRUE)

PCA_tb[PCA_tb[,"ID"] == "alpha54", "Clades"] = "Lab"
PCA_tb[PCA_tb[,"ID"] == "GIL", "Clades"] = "Lab"
PCA_tb[PCA_tb[,"ID"] == "Lang", "Clades"] = "Lab"

## change alpha54 to LCTL1
## change GIL to LCTL2
PCA_tb[PCA_tb[,"ID"] == "alpha54", "ID"] = "LCTL1"
PCA_tb[PCA_tb[,"ID"] == "GIL", "ID"] = "LCTL2"



pc1_axis_str <- paste("PC1 (", round(100*percent_vars[1], 2), "%)", sep="")
pc2_axis_str <- paste("PC2 (", round(100*percent_vars[2], 2), "%)", sep="")

color_vec2 <- c( "#0BC346", #"Wine/Euopean"
                "#D55E00", #African palm wine
                "#AA4499", #Asian islands
                "#332288", #Sake
                "#CAB23C", #Asian fermentation
                "#CC6677", #Brazilian bioehtanol
                "#2F30EE", #Mosaic beer
                "#000000", #control 
                "#999999") #Mosaic region 1



setEPS()
postscript("../plot/fig4_PCA.eps", width=6, height=5.6)

p_PCA <- ggplot(PCA_tb, aes( PC1, PC2, label= ID, color= Clades)) + theme_bw() +
  labs(y= pc2_axis_str, x=pc1_axis_str, size=11) +
  scale_color_manual(values= color_vec2, labels= c("Wine/Euopean",
                                                  "African palm wine", "Asian islands", "Sake", "Asian fermentation",
                                                  "Brazilian bioethanol", "Mosaic beer", "Control", "Mosaic region 1")) +
  scale_shape_manual(values = c(16,5)) +
  theme (panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.title=element_text(size=14,face="bold"),
         axis.text.x = element_text(angle = 45, vjust=0.7, size=12),
         axis.text.y = element_text( size=12),
         legend.text = element_text( size = 11)) + geom_text_repel(size=3, colour = "black")  + geom_point()

p2 <- plot_PCA_vec(mut.pca, c(1,1))

plot_grid(p_PCA, p2,  align = "v", nrow = 2, rel_heights = c(2/3, 1/3))

dev.off()

#p_denovo <- plot_PCA_analysis_tb_raw(mut.pca , "de novo mutation spectrum PCA", c(1,1))
#ggsave("../plot/mut_freq_new_10_21_20_PCA.pdf", width=7, height=5)

##############
# fig 4B. 
##############

## only plot bar plot with three strains 
strains3 <- c("alpha54_316",  "AEQ_252", "AAR_316") 
mut_tb_sub_frq_long_3 <- mut_tb_sub_frq_long[mut_tb_sub_frq_long[,"variable"] %in% strains3 ,] 


p <- ggplot(data=mut_tb_sub_frq_long_3, aes(x=mut, y=value, fill=mut))+ geom_bar(stat="identity") +
  scale_fill_manual( values= color_vec) +facet_wrap(~variable, ncol=4 ) +
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
         axis.text = element_text(size=12, face="bold"),
         legend.position = "none") +
  labs(y= "Frequency")

ggsave("../plot/fig4B_bar_plot.pdf", width=6, height=2.5)


######## 
## plot proportions of MNMs per strain, suppl S7.
########
double_MNM <- read.table("../MNM/final_MNM_double.txt", header=TRUE)
complex_MNM <- read.table("../MNM/final_MNM_complex.txt", header =TRUE)
single <- read.table("../all_mut_wd_excl_MNM/all_mut_combined.txt", header=TRUE)

## calculate # of mutants per category per strain
double_MNM_count <- double_MNM %>% group_by(strain) %>% summarise(n()) %>% data.frame
complex_MNM_count <- complex_MNM %>% group_by(strain) %>% summarise(n()) %>% data.frame
single_count <- single %>% group_by(strain) %>% summarise(n()) %>% data.frame

names(double_MNM_count) <- c("strain", "double_MNM")
names(complex_MNM_count) <- c("strain", "complex_MNM")
names(single_count) <- c("strain", "single")

all_mut <- merge(double_MNM_count, complex_MNM_count, by="strain", all=TRUE)
all_mut <- merge(all_mut, single_count, by="strain", all=TRUE)

all_mut[is.na(all_mut)] <- 0
rownames(all_mut) <- all_mut[,1]
all_mut[,1] <- NULL
all_mut_t <- t(all_mut)

colname <- colnames(all_mut_t)
colname [colname=="alpha54"] <- "LCTL1"
colname [colname=="GIL"] <- "LCTL2"
colnames(all_mut_t) <- colname

all_mut_frq <- count_2_freq_multi_col_all2 (all_mut_t)
## remove frq for single mutations 
all_mut_frq2 <- all_mut_frq[ !rownames(all_mut_frq) %in% c("single"), ]

## get the strains in the order that the MNM mut frq from high to low
strain_order <- names(sort(-colSums(all_mut_frq2)))

## then plot stacked bar of mutation types 
all_mut_frq_long <- reshape2::melt(all_mut_frq2)
names(all_mut_frq_long) <- c("mut", "strain", "frq")
all_mut_frq_long[,"strain"] <- factor(all_mut_frq_long[,"strain"], levels = strain_order)
all_mut_frq_long[,"mut"] <- factor(all_mut_frq_long[,"mut"], levels = c("complex_MNM", "double_MNM"))


color_vec3 <- c("#E69F00", "#0072B2")

p <- ggplot(data=all_mut_frq_long, aes(x=strain, y=frq, fill=mut))+ geom_bar(position= "stack", stat="identity") +
  scale_fill_manual( values= color_vec3) +  
  #facet_grid(~ category, scales = "free_x", space = "free_x") + 
  geom_col(alpha = 0.8, width = 0.85) +
  theme_bw() +
  theme( strip.placement = "outside",
         strip.background = element_blank(),
         strip.text = element_text( size = 11 ),
         panel.border=element_blank(),
         axis.line = element_line(),
         panel.spacing = unit(2, "lines"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 90),
         axis.title=element_text(size=11,face="bold"),
         axis.text = element_text(size=9, face="bold")) +
  labs(y= "Frequency", x="", fill="Mutation") 


ggsave("../plot/figS_mut_frq_per_strain.pdf", width=4.3, height=2.7)



###################
### fig 4A
### plot scatter plot of single mut# rate versus MNM# rate, use mean mutation rate estimate 
###################

## read in mean CAN1 rate for each strain 
## then normalize MNM frq with rate (to get the estimate count)
can1_rate <- read.table("../data/mean_mut_rate_CAN1_strains.txt", header=TRUE)


MNM_frq <- all_mut_frq["double_MNM",] + all_mut_frq["complex_MNM",]
all_mut_frq <- rbind(MNM_frq, all_mut_frq["single",])
rownames(all_mut_frq) <- c("MNM", "single")
all_mut_frq_t <- t(all_mut_frq)
all_mut_frq_df <- data.frame(all_mut_frq_t)
all_mut_frq_df[,"ID"] <- rownames(all_mut_frq_df)

## add annotation 
all_mut_frq_df <- merge(all_mut_frq_df, PCA_tb, by="ID")

## adjust to normalized rate 
names(can1_rate)[1] <- "ID"
all_mut_frq_df_new <- merge(all_mut_frq_df, can1_rate, by="ID")
      

all_mut_frq_df_new[,"MNM_rate"] <- all_mut_frq_df_new[,"MNM"] *all_mut_frq_df_new[,"mean_mut_rate"] 

all_mut_frq_df_new[,"single_rate"] <- all_mut_frq_df_new[,"single"] *all_mut_frq_df_new[,"mean_mut_rate"] 


p2 <- ggplot(data=all_mut_frq_df_new, aes(x=single_rate, y=MNM_rate, label= ID, color= Clades ))+ geom_point() +
  geom_smooth(method=lm,aes(group=1), fullrange= TRUE,colour="black") +
  scale_color_manual(values= color_vec2, labels= c("Wine/Euopean",
                                                   "African palm wine", "Asian islands", "Sake", "Asian fermentation",
                                                   "Brazilian bioethanol", "Mosaic beer", "Control", "Mosaic region 1")) +
  #scale_fill_manual( values= color_vec3) +  
  #facet_grid(~ category, scales = "free_x", space = "free_x") + 
  #geom_col(alpha = 0.8, width = 0.85) +
  theme_bw() +
  theme( 
         panel.border=element_blank(),
         axis.line = element_line(),
         panel.spacing = unit(2, "lines"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.title=element_text(size=11,face="bold"),
         axis.text = element_text(size=9, face="bold")) +
  labs(y="MNM rate", x="Single-bp mutation rate") + geom_text_repel(size=3, colour = "black")  

ggsave("../plot/Fig_4A_MNM_single_scatter.pdf", width=5, height=3)


all_mut[,"MNM"] <- all_mut[,"double_MNM"] + all_mut[,"complex_MNM"]
all_mut[,"ID"] <- rownames(all_mut)

## perform chi-sq test agaist LCTL1, of the proportion of MNM to single mutation
## p<0.05: AAB, CEV, LCTL2. 
for(i in 1:nrow(all_mut))
{
  tb <- rbind(x= all_mut[i,c("MNM", "single")], y= all_mut["alpha54",c("MNM", "single")])
  print( all_mut[i,"ID"])
  print(chisq.test(tb))
}

