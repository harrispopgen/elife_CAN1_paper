

source("est_m_basic.R")

library(scales)
library(ggpubr)

mut_tb <- calc_mut_rate_all("fluc_mut_count_summary_CAN1.csv", "fluc_cell_count_all.csv")

# exclude  CEI-d1-b1

mut_tb_new <- mut_tb [ (mut_tb[,"Strain"] != "CEI-d1-b1"),]

## change alpha54 to LCTL1 
mut_tb_new[ mut_tb_new[,"Strain"]=="alpha54" ,"Strain"] = "LCTL1"
mut_tb_new[ mut_tb_new[,"Strain"]=="GIL" ,"Strain"] = "LCTL2"

## color boxplot with the clade 

annot_tb <- read.table("strain_anno_group.txt", header=TRUE, sep=",",
                       stringsAsFactors = FALSE)

names(annot_tb)[1] <- "Strain"

mut_tb_new <- merge(mut_tb_new,annot_tb, by="Strain", all.x=TRUE)

mut_tb_new [mut_tb_new[,"Strain"] == "LCTL1", "Clades"] = "Control"
mut_tb_new [mut_tb_new[,"Strain"] == "LCTL2", "Clades"] = "Control"

mut_tb_new["isControl"] <- FALSE
mut_tb_new[mut_tb_new[,"Clades"] == "Control","isControl"] <- TRUE

color_vec <- c( "#0BC346", #"Wine/Euopean"
                "#D55E00", #African palm wine
                "#AA4499", #Asian islands
                "#332288", #Sake
                "#CAB23C", #Asian fermentation
                "#CC6677", #Brazilian bioehtanol
                "#2F30EE", #Mosaic beer
                "#000000", #control 
                "#999999") #Mosaic region 1



ylabel <- expression(paste(bold("Mutation rate of "),  bolditalic("CAN1"), bold(" (log 2)")))

p <- plot_mut_rate(mut_tb_new, "Clades", "isControl") + theme_bw() +
  labs(y= ylabel, x="Strain", size=11) +
  scale_color_manual(values= color_vec, labels= c("Wine/Euopean",
                                                  "African palm wine", "Asian islands", "Sake", "Asian fermentation",
                                                  "Brazilian bioethanol", "Mosaic beer", "Control", "Mosaic region 1")) +
  scale_shape_manual(values = c(16,5)) +
  #scale_y_continuous(breaks=c(2e-7, 5e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-6)) +
  #scale_y_continuous(trans = 'log2', labels = scales::scientific) + # log 10 scale ; https://stackoverflow.com/questions/42323247/how-to-force-axis-values-to-scientific-notation-in-ggplot
  scale_y_continuous(trans = "log2", n.breaks = 10,labels = label_scientific(digits = 2))+
  theme (axis.title=element_text(size=14,face="bold"),
         axis.text.x = element_text(angle = 45, vjust=0.7, size=12),
         axis.text.y = element_text( size=12),
         legend.text = element_text( size = 11)) + annotation_logticks(sides="l") 


#ggsave("fig3_by_rate.pdf", width=8.2, height=4.2)
ggsave("fig3_raw.svg", width=8.4, height=4.2)


## calculate mean mutation rate for each strain

mut_mean <- mut_tb_new %>% group_by(Strain) %>% summarise(mean(m_rate)) %>% data.frame 

names(mut_mean) <- c("Strain", "mean_mut_rate")

write.table(mut_mean, "mean_mut_rate_CAN1_strains.txt",sep="\t", quote=FALSE, row.names=FALSE)


## plot mean mutation rate against growth rate measurement in YPD 30C

growth_rate <- read.table("Phenotyping_data_size.txt", header =TRUE, sep="\t")

growth_rate_cols <- growth_rate[,c("Standardized_name", "YPD_40h")]
names(growth_rate_cols) <- c("Strain", "growth")

mut_rate_growth_tb <- merge(mut_mean, growth_rate_cols, by="Strain")

#p_growth <- ggplot(mut_rate_growth_tb, aes( growth, mean_mut_rate)) + geom_point()


## calculate r^2
#r2 <- summary(lm(mut_rate_growth_tb$mean_mut_rate ~ mut_rate_growth_tb$growth ))$r.squared

#r2_label <- parse(text= paste("R^2=", round(r2, digits=5), sep=""))
#p_growth <- p_growth+ geom_text(x=1000, y=2e-6, label=r2_label)

ggscatter(mut_rate_growth_tb, x = "growth", y = "mean_mut_rate") +
  stat_cor(label.y = 2e-6,
           label.x = 950,
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  expand_limits(x = 1200) +
  labs(y= "mean mutation rate", x="growth rate (colony size)", size=11)

ggsave("mut_rate_growth.pdf", width=5.5, height=4.5)


