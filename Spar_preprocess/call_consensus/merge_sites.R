library(reshape2)

## modified on 4.29.19, change the output into bed file
## merge SNPs from different S.para, only include 
merge_sites_Spara <- function()
{
	S_para_list = c( "YPS138", "UWOPS919171", "N44", "UFRJ50816")
	
	file_name <- paste("data/CBS432", "sites", "excl_indel.bed", sep="_")
	tb<-read.table(file_name)
	tb[,"key"] <- paste(tb[,1],tb[,3],sep="_") ## merge chr and postion (use 1-based)
	merged_tb <- tb [,c("V4","key")]
	names(merged_tb) <- c("CBS432","key")
	
	for (S_para in S_para_list)
	{
		file_name <- paste("data/", S_para, "_sites_excl_indel.bed", sep="")
		tb<-read.table(file_name)
		tb[,"key"] <- paste(tb[,1],tb[,3],sep="_")
		subtb <- tb [,c("V4","key")]
		names(subtb) <- c(S_para,"key")
		merged_tb <- merge(merged_tb, subtb, by ="key", all=TRUE)
	}
	
	## only keep SNPs that have at most 1 NA and the others having the same ALT allele
	at_most_1_miss <- merged_tb [rowSums(is.na(merged_tb)) <=1 ,]
	
	# count the unique number of values for each row 
	# https://stackoverflow.com/questions/2269084/handling-na-values-in-apply-and-unique
	count_alt_allele <- apply(at_most_1_miss, 1, function(x)length(unique(x[!is.na(x)]))) ## count for each line, unique values (key also counts as one value)
	shared_alt_allele <- at_most_1_miss[count_alt_allele==2, ] ## count=2, means there are only one type of value
	## output all the sites that are at most 1 missing value, and all the outgroup have the same allele at the site
	write.table(shared_alt_allele, "data/merged_shared_alt_SNP_para.txt",sep="\t",quote=FALSE, row.names=F)
	
	# extract unique values (alt alleles) for each site
	alt_allele <- apply(shared_alt_allele[,c(2,3,4,5,6)], 1, function(x) unique(x[!is.na(x)]))
	
	alt_allele_tb <- colsplit(shared_alt_allele[,1],pattern="_", names= c("chr","pos"))
	alt_allele_tb[,"begin"] <- alt_allele_tb[,"pos"] -1 
	alt_allele_tb[,"allele"] <- alt_allele
	write.table(alt_allele_tb[,c("chr","begin","pos","allele")], "data/Scer_para_consol_SNPs.txt",sep="\t",quote=FALSE, row.names=F, col.names=F) ## no headers
	
}

merge_sites_Spara()
