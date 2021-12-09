## plot phylogenetic tree of y1002 
## and color the nodes with p<0.05 in a different color 

library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(colorspace)
library(tibble)
library(tidytree)


## change this script to take arguments from command line
## ==1: use the previous 5% cutoff (C>A 5% outlier)
## ==2: use the 5% overlap with the bootstraping method 




args = commandArgs(trailingOnly=TRUE)
if(length(args)==1)
{
  flag= strtoi(args[1])
} else
{
  print ("error! please specify argument")
  quit(save="no")
}


## use phylogenetic tree built from SNPs 
y1002_dist <- read.table("/net/harris/vol1/project/yeast_1002/nobackup/1011DistanceMatrixBasedOnSNPs.tab",header=TRUE,row.names=1)


## generate a subset of tree using strains present in the file (excluding closely related strains)
CA_ratio_pval_tb <- read.table("../data/y1002_C_A_ratio_pval.txt", header=TRUE, stringsAsFactors = FALSE)
CA_ratio_pval_tb_sub <- CA_ratio_pval_tb[,c("strains", "pval")]
names(CA_ratio_pval_tb_sub) <- c("strains", "pval1")

## read in the second file with simulations 
CA_ratio_pval_tb_simu2 <- read.table("../data/y1002_C_A_ratio_pval_newsimu.txt", header=TRUE, stringsAsFactors = FALSE)
CA_ratio_pval_tb_simu2_sub <- CA_ratio_pval_tb_simu2[,c("strains", "pval")]
names(CA_ratio_pval_tb_simu2_sub)  <- c("strains", "pval2")

CA_ratio_pval_tb_merge <- merge(CA_ratio_pval_tb_sub, CA_ratio_pval_tb_simu2_sub, by="strains")

subset <- CA_ratio_pval_tb_merge[,"strains"]


## exclude strain "ARC", not in the distance matrix 
idx2 <- subset!="ARC"
subset2 <- subset[idx2]


CA_ratio_pval_tb2 <- CA_ratio_pval_tb_merge[idx2,]

## will color subset of strains with pval <0.05
idxA <- CA_ratio_pval_tb2[,"pval1"] <0.05
idxB <- CA_ratio_pval_tb2[,"pval2"] <0.05
if(flag==1)
{
  idx <- idxA
}else if (flag==2)
{
  idx <- idxA & idxB
}



CA_ratio_pval_tb2[,"group"] <- NA
CA_ratio_pval_tb2[idx,"group"] <- "1"

## write group file as the annotation file to plot the tree 
subtb <- CA_ratio_pval_tb2[,c("strains", "group")]



## read in annotation of strains 
names(subtb) <- c("ID", "CAgroup")
annot <- read.table("/net/harris/vol1/project/yeast_1002/data/strain_anno_group.txt", header =TRUE, sep=",",
                    stringsAsFactors = FALSE)
annot_new <- annot[annot[,"group"] != "-1",]


## strains that has group information 
hasgroup_strain <- annot_new[,"ID"]
subset3 <- subset2[subset2 %in% hasgroup_strain]

y1002_dist_sub <- y1002_dist[subset3,subset3]


subtb_annot <- merge(subtb, annot_new, by="ID")



## then generate a NJ tree by this subset of strains
y1002_sub_nj <- NJ(as.matrix(y1002_dist_sub))

## add pseudo root
y1002_sub_nj_new <- root(y1002_sub_nj, "CEG")



## only keep Clades and group
subtb_annot2 <- subtb_annot[,c("ID", "CAgroup", "Clades", "group")]


names(subtb_annot2) <- c("label", "CAgroup", "Clades", "group")

subtb_annot2new <- subtb_annot2[,c("label", "CAgroup")]



## then add group information by clade 

## exclude strains which are not assignned to a clade
subtb_annot3 <- subtb_annot2[!is.na(subtb_annot2[,"Clades"]),]

## trim spaces 
subtb_annot3[,"Clades"] <- trimws(subtb_annot3[,"Clades"])

names(subtb_annot3) <- c("label", "CAgroup", "Clades", "group")


## https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ggtree/inst/doc/treeManipulation.html
## generate a list of tip nodes, to feed into ggtree, to color sub-branches
clades <- levels(factor(subtb_annot3[,"Clades"]))
nclades <- length(clades)
clades_for_ggtree <- list()
for(i in 1:nclades)
{
  clade_name <- clades[i]
  subset <- subtb_annot3 [subtb_annot3[,"Clades"] == clade_name,"label"]
  clades_for_ggtree[[clade_name]] <- subset
}



tree <- groupOTU(as_tibble(y1002_sub_nj_new), clades_for_ggtree)




## add annotation of CA group
tree2 <- as.treedata(full_join(tree, subtb_annot2new, by="label"))


## find the node for each clade 
## 1. wine 
n1 <- MRCA(tree2, c("ATM", "AHD"))
## 3. brazilian ethanol
n2 <- MRCA(tree2, c("BVC", "CNL"))
## 4. medteranina oak
n3 <- MRCA(tree2, c("BFP", "CCL"))
## 5. french dairy 
n4 <- MRCA(tree2, c("ARS", "BGE"))
## 6. african beer 
n5 <- MRCA(tree2, c("ANM", "AFL"))
## 9.mexican agave
n6 <- MRCA(tree2, c("CPL", "CPN"))
## 10. french guaina human 
n7 <- MRCA(tree2, c("BNB", "BDQ"))
## 11. ale beer 
n8 <- MRCA(tree2, c("CGC", "CFC"))
## 12. west african coca
n9 <- MRCA(tree2, c("CQN", "CQD"))
## 13. african palm wine 
n10 <- MRCA(tree2, c("AKH", "BAA"))

## 24. asian islands 
n11 <- MRCA(tree2, c("CCV", "ARH"))
## 25. sake
n12 <- MRCA(tree2, c("ARP", "ANR"))
## 26. asian ferm 
n15 <- MRCA(tree2, c("AKL", "AHQ"))

##n13 <- MRCA(y1002_sub_nj_new, c("AST", "BHL"))

##n14 <- MRCA(y1002_sub_nj_new, c("CDM", "AER"))




## zoom in to scale the clade that has specific strains 
z1 <- MRCA(tree2, c("SACE_YAG", "BRM"))
z2 <- MRCA(tree2, c("CHG", "AFF"))
z3 <- MRCA(tree2, c("ADT", "ACN"))

## 7. mosaic beer 
n17 <- MRCA(tree2, c("SACE_YAG", "CBN"))
## 8. mixed origin 
n18 <- MRCA(tree2, c("CHG", "AFF"))
## M2. mosaic region 2
n16 <- MRCA(tree2, c("BDD", "BDK"))




## 
treetb <- as_tibble(tree2)
n14_pre <- treetb[treetb[,"label"] == "CEG","node"]
n14 <- n14_pre[!is.na(n14_pre)]

color_vec <- c( "1. Wine/European"= "#5B1A18", 
                "1. Wine/European (subclade 1)"= "#5B1A18",
                "1. Wine/European (subclade 2)"= "#5B1A18",
                "1. Wine/European (subclade 4)" = "#5B1A18",
                "2. Alpechin" = "#f76a8c",
                "3. Brazilian bioethanol"="#f76a8c",
                "4. Mediterranean oak" = "#377eb8",
                "5. French dairy"= "#fcae91", 
                "6. African beer"=  "#377eb8", 
                "7. Mosaic beer"="#6ADD77", 
                "8. Mixed origin"= "#4C00A3",
                "9. Mexican agave" = "#6ADD77",
                "10. French Guiana human" = "#fcae91",
                "11. Ale beer" = "#FAD510",
                "12. West African cocoa" = "#FAD510",  
                "13. African palm wine"= "#9c5518",
                "14. CHNIII" = "#CC6677",
                "15. CHNII" = "#CC6677",
                "16. CHNI" = "#CC6677",
                "17. Taiwanese" =  "#882255",
                "18. Far East Asia" = "#CC6677",
                "19. Malaysian" = "#CC6677",
                "20. CHN V" = "#CC6677",
                "21. Ecuadorean" = "#CC6677",
                "22. Far East Russian" = "#CC6677",
                "23. North American oak"="#3CAEA3",
                "24. Asian islands"= "#fbcffc",
                "25. Sake" = "#be79df", 
                "26. Asian fermentation"= "#2b580c",
                "M1. Mosaic region 1" = "#be79df",
                "M2. Mosaic region 2" = "#f76a8c",
                "M3. Mosaic region 3" = "#87bbd6",
                "0" = "#999999",
                "1" = "#C70405"
)


## drop certain tips 

#drop.tip(as.treedata(tree2), to_drop)


p <- ggtree(tree2, layout= "circular", size=0.6) 
   
p <- rotate_tree(p, 180)

p <- p + geom_tiplab(aes(color= CAgroup),  size=8, hjust = 1)  +
  scale_color_manual(values=color_vec)

#p <- p + geom_tiplab(  size=3, hjust = 1)  +
#  scale_color_manual(values=color_vec)

p <- scaleClade(p, n1, 0.3)
p <- scaleClade(p, n2, 0.5)
p <- scaleClade(p, n3, 0.65)
p <- scaleClade(p, n4, 0.5)
p <- scaleClade(p, n5, 0.5)
p <- scaleClade(p, n6, 0.5)
p <- scaleClade(p, n7, 0.5)
p <- scaleClade(p, n8, 0.5)
p <- scaleClade(p, n9, 0.5)
p <- scaleClade(p, n10, 0.5)
p <- scaleClade(p, n11, 0.5)
p <- scaleClade(p, n12, 0.5)
#p <- scaleClade(p, n14, 0.2)
p <- scaleClade(p, n15, 0.5)
#p <- scaleClade(p, n16, 0.55)

p <- scaleClade(p, z1, 5.2)
p <- scaleClade(p, z2, 5.2)
p <- scaleClade(p, z3, 5.2)

#p2 <- p %>% collapse(node=n1, mode= "mixed", fill="darkgreen") # +
# # geom_point2(aes(subset=(node==n1)), shape=21, size=3, fill='green')


#p2 <- p %>% scaleClade(n1, 0.001) %>% collapse(n1, 'mixed', fill="darkgreen")  

#p2 <-  collapse(p2, node=n2, mode= "mixed", fill="darkgreen") # +
#  geom_point2(aes(subset=(node==n2)), shape=21, size=3, fill='green')

#p2 <- scaleClade(p, n2, 0.2) %>% collapse(n2, 'min', fill="blue")  

p2 <- collapse(p, n1, 'min', fill="#999999")  

p2 <- collapse(p2, n2, 'min', fill="#999999")  

p2 <- collapse(p2, n3, 'min', fill="#999999")  

p2 <-  collapse(p2, n4, 'min', fill="#999999") 

p2 <- collapse(p2, n5, 'min', fill="#999999")  

p2 <- collapse(p2, n6, 'min', fill="#999999")  

p2 <- collapse(p2, n7, 'min', fill="#999999")  

p2 <-  collapse(p2, n8, 'min', fill="#999999") 

p2 <- collapse(p2, n9, 'min', fill="#999999")  

p2 <- collapse(p2, n10, 'min', fill="#999999")  

p2 <-  collapse(p2, n11, 'min', fill="#999999") 

p2 <- collapse(p2, n12, 'min', fill="#999999")  


p2 <-  collapse(p2, n15, 'min', fill="#999999") 
#p2 <-  collapse(p2, n16, 'min', fill="blue") 


#p2 <-  collapse(p2, node=n11)  +
#  geom_point2(aes(subset=(node==n11)), shape=21, size=3, fill='green')


p2 <- p2 + geom_cladelabel(n1, "Wine/European", fontsize=8, offset=0, hjust=1.2) +
  geom_cladelabel(n2, "Brazilian bioethanol", angle=90, fontsize=8, offset=0.7, hjust=0.3) +
  geom_cladelabel(n3, "Mediterranean oak",  angle=90, fontsize=8, offset=0.5, hjust=0) +
  geom_cladelabel(n4, "French dairy", angle=90, fontsize=8, offset=0.8, hjust=0.5) +
  geom_cladelabel(n5, "African beer", angle=90, fontsize=8, offset=0.6, hjust=0.3) +
  geom_cladelabel(n17, "Mosaic beer",angle=10, fontsize=8, offset=0.5, hjust=1, barsize=1) + 
  geom_cladelabel(n18, "Mixed origin",angle=30, fontsize=8, offset=0.5, hjust=1, barsize=1) + 
  geom_cladelabel(n6, "Mexican agave",angle=300, fontsize=8, offset=0.8, hjust=0.5) + 
  geom_cladelabel(n7, "French Guiana human",angle=300, fontsize=8, offset=0.8, hjust=0.5) + 
  geom_cladelabel(n8, "Ale beer",angle=300, fontsize=8, offset=0.8, hjust=1) + 
  geom_cladelabel(n9, "West African cocoa",angle=270, fontsize=8, offset=0.3, hjust=0) + 
  geom_cladelabel(n10, "African palm wine",angle=10, fontsize=8, offset=0.9, hjust=0.5) + 
  geom_cladelabel(n11, "Asian islands",angle=10, fontsize=8, offset=0.9, hjust=0.5) + 
  geom_cladelabel(n12, "Sake",angle=40, fontsize=8, offset=0.2, hjust=1) + 
  geom_cladelabel(n15, "Asian fermentation",angle=40, fontsize=8, offset=0.2, hjust=1) + 
  geom_cladelabel(n14, "Taiwanese",angle=0, fontsize=8, offset=0.1, hjust=1) + 
  #geom_cladelabel(n16, "Mosaic region 2",angle=0, fontsize=8, offset=0.5, hjust=1) + 
  geom_hilight(node=n17, fill="#2F30EE", alpha=0.45) +
  geom_hilight(node=n18, fill="#FFEE08", alpha=0.45) +
  geom_hilight(node=n14, fill="#999999", alpha=0.5) +
  geom_strip("CIH","ANA", label = "Mosaic region 3", fontsize=8,barsize=1, hjust=1) +
  geom_strip("SACE_GAV","CNC",barsize=1) 

if(flag==1)
{
  ggsave("../plot/fig5_tree.svg", width=21,height=21)
}else if (flag==2)
{
  ggsave("../plot/figS_tree_v2.pdf", width=21,height=21)
}



################
### alternative method using groupOTU




## https://aschuerch.github.io/posts/2017-04-24-blog-post-1
## https://groups.google.com/forum/#!topic/bioc-ggtree/o35PV3iHO-0
#p %<+% CA_ratio_pval_tb2 +
#  geom_tiplab2( aes(color= factor(CAgroup)),
#                size=3)  


# geom_tiplab(CA_ratio_pval_tb2, aes(fill= factor(group)),
#              color = "black",
#              geom = "label")
