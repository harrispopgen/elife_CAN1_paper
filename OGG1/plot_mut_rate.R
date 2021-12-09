
setwd("/Users/Pengyao/Google Drive/yeast_mutator/fluctuation")

source("est_m_basic.R")

library(scales)
library(ggpubr)
library(wesanderson)


mut_tb <- calc_mut_rate_all("can1_3_19_ogg1_plasmids_mut.csv", "ogg1_plasmids_cell_counts.csv")


## create a df to annotate strains
col1 <- c("20510_h2c_emp", "20510_h2c_AEQ_A", "20510_h2c_AEQ_B", "20510_h2c_alpha54_A", "20510_h2c_alpha54_B")
OGG1_AEQ_A_name <- expression(italic("OGG1")^"AEQ/AAR"*"_A")
OGG1_AEQ_B_name <- expression(italic("OGG1")^"AEQ/AAR"*"_B")
OGG1_wt_A_name <- expression(italic("OGG1")^"wt"*"_A")
OGG1_wt_B_name <- expression(italic("OGG1")^"wt"*"_B")
col2 <- c("empty_plasmid", "OGG1_AEQ_A", "OGG1_AEQ_B", "OGG1_wt_A", "OGG1_wt_B")
col2_display <- c("empty_plasmid", OGG1_AEQ_A_name, OGG1_AEQ_B_name, OGG1_wt_A_name, OGG1_wt_B_name)
col3 <- c("empty plasmid", "OGG1_AEQ", "OGG1_AEQ", "OGG1_wt","OGG1_wt")
OGG1_AEQ_name <- expression(italic("OGG1")^"AEQ/AAR")
OGG1_wt_name <- expression(italic("OGG1")^"wt")
col3_display <- c("empty plasmid", OGG1_AEQ_name, OGG1_AEQ_name, OGG1_wt_name, OGG1_wt_name)


df <- data.frame (Strain = col1, Strain_name=col2, Type=col3 )
mut_tb_new <-merge(mut_tb, df, by ="Strain")


ylabel <- expression(paste(bold("Mutation rate of "),  bolditalic("CAN1")))
mut_tb_new[,"Strain"] <-  mut_tb_new[,"Strain_name"] 
mut_tb_new[,"Strain"] <- factor(mut_tb_new[,"Strain"], levels= col2)



p <- ggplot (mut_tb_new, aes_string("Strain","m_rate")) + geom_boxplot(aes_string(color= "Type")) +
  geom_point(aes_string(color= "Type"), size=1.8) +
  theme_bw() +
  scale_color_brewer(palette="Dark2", breaks= col3, labels=col3_display, name="Plasmid") +
  labs(y= ylabel, x="Strain", size=11) +
  scale_y_continuous(trans = "log2", n.breaks = 10,labels = label_scientific(digits = 2))+
  theme (axis.title=element_text(size=14,face="bold"),
         axis.text.x = element_text(angle = 30, vjust=0.7, size=12),
         axis.text.y = element_text( size=12),
         legend.text = element_text( size = 11),
         legend.text.align = 0) + annotation_logticks(sides="l") +
  scale_x_discrete(labels=col2_display) #+ ## change x axis display name 
  #scale_fill_discrete(breaks= col3, labels=col3_display) ## change legend labels


#ggsave("fig6_mut_rate_OGG1_plasmid.svg", width=6, height=4)
#ggsave("mut_rate_OGG1_plasmid.pdf", width=6, height=4)
ggsave(p, file="fig6_mut_rate_OGG1_plasmid.eps", device="eps" , width=6, height=4)



## calculate mean mutation rate for each strain

mut_mean <- mut_tb %>% group_by(Strain) %>% summarise(mean(m_rate)) %>% data.frame 
mut_count_mean <- mut_tb %>% group_by(Strain) %>% summarise(mean(m)) %>% data.frame 


names(mut_mean) <- c("Strain", "mean_mut_rate")

write.table(mut_mean, "mean_mut_rate_CAN1_strains.txt",sep="\t", quote=FALSE, row.names=FALSE)


