## count the overall percentage of different kinds of mutations 

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)




tb <- read.table("../data/all_mut_combined_coord_converted.txt", header=TRUE)


## combine tb with annotation (use full tb)

## read in the annotation file (output from VEP) 
## three types: missense_variant; stop_gained; synonymous_variant
mut_annot <- read.table("../data/mut_VEP_CAN1.csv", sep=",", header =TRUE)
mut_annot_subcol <- mut_annot[,c("Uploaded_variation", "Consequence")]

tb[,"mut"] <- paste(tb[,"ref"], tb[,"variable"], sep="/")
tb[,"Uploaded_variation"] <- paste("V", tb[,"genome_pos"], tb[,"mut"],sep="_" )

## annotate the other non-ref site to be "special"
tb[,"mut_and_pos"] <- paste(tb[,"new_pos"], tb[,"ref"], tb[,"variable"], sep="_" )
tb[tb[,"mut_and_pos"] == "463_G_A", "special"] = "yes"

## only annotate those that have reference 
tb_mut_all_annot <- merge(tb, mut_annot_subcol, by="Uploaded_variation", all.x=TRUE)


## write to file 
write.table(tb_mut_all_annot, "../data/all_mut_combined_coord_annot.csv", quote=FALSE, sep="," , row.names=FALSE)


## get only SNV mutations for annotation
single_mut <- tb[(tb[,"variable"] != "insertions") & (tb[,"variable"] != "deletions") , ]

## exclude one mutation which occurs at the insertion seq, not in ref. 
single_mut2 <- single_mut[single_mut[,"special"] !="yes",]

single_mut2[,"mut_and_pos"] <- paste(single_mut2[,"new_pos"], single_mut2[,"ref"], single_mut2[,"variable"], sep="_" )

## exclude one mutation in AEQ which has a different ref 
single_mut2_excl <- single_mut2[single_mut2[,"mut_and_pos"] != "463_G_A",]

single_mut2_excl[,"mut"] <- paste(single_mut2_excl[,"ref"], single_mut2_excl[,"variable"], sep="/")

single_mut2_excl[,"Uploaded_variation"] <- paste("V", single_mut2_excl[,"genome_pos"], single_mut2_excl[,"mut"],sep="_" )


single_mut2_annot <- merge(single_mut2_excl, mut_annot_subcol, by="Uploaded_variation")

########
#### new: check the occurrence of mutations per site 
## calculate how many of the mutants occur at a unique sites
count_by_pos <- single_mut2_annot %>% group_by(genome_pos)  %>% summarise(sum(count)) %>% data.frame
names(count_by_pos) <- c("genome_pos", "Count")
unique_mut <- count_by_pos[count_by_pos[,2] ==1,]
## how many are unique sites 
nrow(unique_mut)
brk <-seq(0,105)
mut_pos_count_hist <- hist(count_by_pos[,2], breaks=brk, right=TRUE)
## location: 32399 has >100 mutations per site
count_by_pos[count_by_pos[,2] >100,]

## count by different types of mutations at 32399
hotspot1 <- single_mut2_annot[single_mut2_annot[,"genome_pos"] == 32399,]
## G/C: 19, G/T: 85
hotspot1 %>% group_by(mut)   %>% summarise(sum(count))


## output  the count by pos 
write.table(count_by_pos, "../data/CAN1_mut_count_by_occurrence.txt",sep="\t", quote=FALSE,row.names=FALSE)


## will just plot count by position 
geom_pos <- count_by_pos[,"genome_pos"]
p_count_by_pos <-  ggplot(count_by_pos,aes(x=genome_pos, y=Count)) + geom_point( color= "#69b3a2") +
  geom_linerange(aes(x=genome_pos, ymax=Count, ymin=0), color= "#69b3a2") +
  geom_text_repel(data=subset(count_by_pos, Count>50), aes(label=genome_pos)) + 
  theme_bw() +
  theme (axis.title=element_text(size=12,face="bold"),
         axis.text.x = element_text(angle = 45, vjust=0.7, size=10),
         axis.text.y = element_text( size=10),
         legend.text = element_text( size = 9))  

ggsave("../plot/FigS_mut_count_by_site.pdf", width=5, height=5)



## calculate binomial distribution (for expected # of mutants per site)
## use total observed sites that get the mutations 

hist_p <- ggplot(count_by_pos,aes(x=Count)) + geom_histogram( binwidth=1, fill="#69b3a2", color="white") +
  theme_bw() +
  theme (axis.title=element_text(size=12,face="bold"),
         axis.text.x = element_text(angle = 45, vjust=0.7, size=10),
         axis.text.y = element_text( size=10),
         legend.text = element_text( size = 9)) +
  xlab("mut count")  + ylab("Count")


ggsave("../plot/FigS_mut_count_per_site_obs.pdf", width=5, height=4)

p_mut <- 1/476
x <- seq(1,105)
y_exp <- dbinom (x, 4561, p_mut) * 4561
exp <- data.frame(x=x, y=y_exp)

hist_p <- ggplot(count_by_pos,aes(x=Count)) + geom_histogram( binwidth=1, fill="#69b3a2", color="white") +
  geom_line (data= exp, aes(x=x,y=y), color= '#7c5295') +
  theme_bw() +
  theme (axis.title=element_text(size=14,face="bold"),
         axis.text.x = element_text(angle = 45, vjust=0.7, size=12),
         axis.text.y = element_text( size=12),
         legend.text = element_text( size = 11))  

ggsave("../plot/mut_count_per_site_add_exp.pdf")

## calculate mutation counts per unique mutation 
count_by_uniq <- single_mut2_annot %>% group_by(Uploaded_variation)  %>% summarise(sum(count)) %>% data.frame
brk2 <-seq(0,85)
mut_count_hist <- hist(count_by_uniq[,2], breaks=brk2, right=TRUE)


## evaluate occurrence in syn mutations
single_mut2_annot_syn <- single_mut2_annot [single_mut2_annot[,"Consequence"] == "synonymous_variant",]
single_mut2_annot_syn_w_count <-merge(single_mut2_annot_syn,count_by_pos, by="genome_pos" )



#######

## missense_variant: 2676
## stop_gained: 1866
## synonymous_variant: 17
single_mut2_annot %>% group_by(Consequence) %>% summarise(sum(count))

strain_mut_type_count <- single_mut2_annot %>% group_by (strain, Consequence)  %>% summarise(sum(count)) %>% data.frame
strain_mut_missense_count <- strain_mut_type_count [strain_mut_type_count[,"Consequence"] == "missense_variant",]
strain_mut_nonsense_count <- strain_mut_type_count [strain_mut_type_count[,"Consequence"] == "stop_gained",]
names(strain_mut_missense_count) <- c("strain", "Consequence", "missense_count")
names(strain_mut_nonsense_count) <- c("strain", "Consequence", "nonsense_count")

strain_mut_merged <- merge(strain_mut_missense_count,strain_mut_nonsense_count, by="strain" )
strain_mut_merged[,"mis2non_ratio"] <- strain_mut_merged[,"missense_count"] /strain_mut_merged[,"nonsense_count"]

ordered_strains <- strain_mut_merged[order(strain_mut_merged[,"mis2non_ratio"]),"strain"]
strain_mut_merged[,"strain"] <- factor(strain_mut_merged[,"strain"] , levels =ordered_strains )

## scatter plot of #of missense and nonsense mutations per strain 
p <- ggplot(strain_mut_merged, aes(x=nonsense_count, y=missense_count, label=strain)) + geom_point() +
  geom_text_repel(size=3, colour = "black")  
ggsave("../plot/type_of_mut_per_strain.pdf")


p <- ggplot(strain_mut_merged, aes(x=strain, y=mis2non_ratio)) + geom_bar(stat="identity") 
ggsave("../plot/type_of_mut_ratio_per_strain.pdf")

## plot the position of mutations per strain 
## combine muts by location
mut_count_per_pos <- single_mut2_annot %>% group_by(strain, new_pos) %>% summarise(sum(count)) %>% data.frame
names(mut_count_per_pos) <- c("strain", "pos", "count")
p <- ggplot(mut_count_per_pos, aes(x=pos, y= count)) + geom_bar(stat="identity") +
  facet_wrap(~strain, nrow=9 ) 
ggsave("../plot/mut_pos_per_strain.pdf", width=8, height=16)

## check if it correlates with mutation rate
mut_rate_can1 <- read.table("../data/mean_mut_rate_CAN1_strains.txt", header =TRUE)
names(mut_rate_can1) <- c("strain", "mean_mut_rate")
strain_mut_merged2 <- merge(strain_mut_merged, mut_rate_can1, by="strain" )

ggscatter(strain_mut_merged2, x = "mis2non_ratio", y = "mean_mut_rate") +
  stat_cor(label.y = 2e-6,
           label.x = 2.2,
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  expand_limits(x = 3) +
  labs(y= "mean mutation rate", x="missense to nonsense mutation ratio", size=11)

ggsave("../plot/cor_mis_non_ratio_to_mut_rate.pdf", width=5.5, height=4.5)


## will calculate total # of unique sites of mutations per strain 

