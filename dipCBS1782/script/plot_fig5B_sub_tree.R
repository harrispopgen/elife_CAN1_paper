## plot subtree of 1011 of some beer strains 

library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(colorspace)
library(tibble)
library(tidytree)
library(treeio)
library(ggrepel)

## use phylogenetic tree built from SNPs 
y1002_dist <- read.table("/net/harris/vol1/project/yeast_1002/nobackup/1011DistanceMatrixBasedOnSNPs.tab",header=TRUE,row.names=1)


## generate a tree using NJ (newick format)
y1002_nj<-NJ(as.matrix(y1002_dist))


## remove AEA 
y1002_nj_new <- drop.tip(y1002_nj, "AEA")

## only get subset of beer strains
strains <- c("BRM", "SACE_YAG", "AEQ", "AAR", "SACE_YAB", "AQH", "AFP")


## get the node for this clade 
node_show <- MRCA(y1002_nj_new, strains)

#p <- ggtree(y1002_nj_new)

# https://tbradley1013.github.io/2018/06/19/subsetting-phylogenetic-trees/
# https://rdrr.io/github/GuangchuangYu/treeio/man/tree_subset.html
subtree <- tree_subset(y1002_nj_new, node = node_show, levels_back =0 )

## define two clades, one has enrichment of C>A, the other don't
groups <- list(c1= c("BRM", "SACE_YAG", "AEQ", "AAR"),
               c2 = c("SACE_YAB", "AQH", "AFP","AFA", "AFB"))

subtree2 <- groupOTU(subtree, groups)


#p1 <- scaleClade(p, node=node_show, scale=2) +geom_tiplab()

#add circle to the end of tip
#https://phe-bioinformatics.github.io/blog/2017/07/07/ggtree_part2
p <- ggtree(subtree2)+geom_tiplab( size=8) + geom_tippoint(aes(color=group), size=3)

#p2 <- open_tree(p, 180) + geom_text_repel(aes(label=label)) + scale_color_manual(values=c("firebrick","black"))

p2 <- p  + scale_color_manual(values=c("#1c7c54","#999999")) + xlim(0,0.8)
#p3 <- viewClade(p+geom_tiplab(), node=node_show)

ggsave("../plot/fig5_subtree.svg", width=5.5, height=4.5)


