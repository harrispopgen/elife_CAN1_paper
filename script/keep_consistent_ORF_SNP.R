# for SNPs that in more than one ORF, if the strand information is consistent, keep it, otherwise remove this SNP

library(dplyr)

snp_tb <- read.table("../data/gene_gff_in_anc_consise.bed")

names(snp_tb) <- c("chr", "begin", "end", "allele", "strand")

snp_tb[,"key"] <- paste(snp_tb[,1], snp_tb[,2], sep="_")

## duplicates 

duplicates <- unique(snp_tb[duplicated(snp_tb[,"key"]),"key"])

## non-duplicates

unique <- snp_tb[! (snp_tb[,"key"] %in% duplicates),]

## get subset of duplicated SNP where the strand is the same 

snp_tb_duplicates <- snp_tb[snp_tb[,"key"] %in% duplicates,]
snp_tb_duplicates_count <- snp_tb_duplicates %>% group_by(key) %>% summarise(length(unique(strand))) %>% data.frame

snp_same_oriet <- snp_tb_duplicates_count [snp_tb_duplicates_count[,2] ==1, 1]

## then get those ones 
resolve_dup <- unique (snp_tb_duplicates [snp_tb_duplicates[,"key"] %in% snp_same_oriet,])


new_snp_tb <- rbind(unique, resolve_dup)

new_snp_tb_order <- new_snp_tb[ order(new_snp_tb$chr, new_snp_tb$begin),]

## drop the key
new_snp_tb_order[,"key"] <- NULL

write.table(new_snp_tb_order, "../data/gene_gff_in_anc_consise_resolve.bed", quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")
