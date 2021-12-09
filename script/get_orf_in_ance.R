## get orf that in ancestral called sites (input file has all the orf, and repeated 
## the number of times that a site is present in ancestral)

library(dplyr)

gene_overlap_tb <- read.table("../data/gene_coord_in_anc.bed")

names(gene_overlap_tb) <- c("chr", "begin", "end", "gene")

## get the coverage of each gene as well as the length of the gene 
gene_overlap_tb %>%  group_by (gene) %>% summarise(n=n(), begin = min(begin), end = min(end), len= min(end)- min(begin)) %>% data.frame  -> gene_cov_count

## only include genes that has >90% coverage
gene_cov_count[,"ratio"] <- gene_cov_count[,"n"]/gene_cov_count[,"len"]
gene_cov_count_sub <- gene_cov_count [gene_cov_count[,"ratio"] >=0.9,]

## write (only gene names) to file 
write.table(gene_cov_count_sub[,"gene"], "../data/gene_in_anc.txt",  quote=FALSE, row.names=FALSE, col.names=FALSE)
