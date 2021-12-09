## extract indels from the outgoup file; will be used by plink to filter out SNPs that are within the regions


## outfile format: chr, begin, end, gene_name (here just use "Indel")
## plink format info: http://zzz.bwh.harvard.edu/plink/dataman.shtml#exclude

## differences in 0 or 1-based for different format
# http://bedtools.readthedocs.io/en/latest/content/overview.html

## generate bed file with only indels 
indel_bed <- function(S_other)
{

	chrs <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

	## note: need to change chromosome name to numbers! (which is the how it is encoded in vcf files)

	## read the chr size file; to replace the last line for each chr
	chr_sizes <- read.table("/net/harris/vol1/data/Pengyao/nobackup/S288C_ref/sacCer3.chrom.sizes", stringsAsFactors=FALSE)


	out_full_tb <- NULL
	out_full_tb_no <- NULL ## chromosome in numbers
	chr_no =1
	for (chr in chrs)
	{
 	 ## read in the SNP_indel diff table per chr
 	 file_name <- paste("data/Scer_", S_other, ".sorted_",chr,"_diff_tb.txt",sep="")
  	 diff_tb<- read.table(file_name, stringsAsFactors=FALSE)
 	 indel_out_tb <- diff_tb[diff_tb[,3]=="Indel",c(1,4,5,3)]
  	
  	 indel_out_tb_bed <- indel_out_tb[,c(1,2,3)]

 	 ##replace the last line in the table (based on actual size)
  	 last_pos <- chr_sizes[chr_sizes[,1]==chr,2]
  	 indel_out_tb_bed[nrow(indel_out_tb_bed),3] <- last_pos

  	 ## exclude any line which is beyond the end of the chromosome
  	 indel_out_tb_bed_new <- indel_out_tb_bed[indel_out_tb_bed[,2] <= last_pos,]

  	 out_full_tb <- rbind(out_full_tb, indel_out_tb_bed_new)
	}
	
	####### important! The bed file used in bedtools and UCSC is 0-based, and ends with 1-based
	## what my code generates for begin and end are all 1-based
	## for exaple the first base: 
	## chr1   0        1    first_base

	out_full_tb[,2] <- as.numeric(out_full_tb[,2]) -1
	
	## replace chrI ..  chrXI to 1..16 and output in bed format (to extract subset of SNPs using bedtools)
	chr_no=1
	new_tb <- NULL
	for (chr in chrs)
	{
 	 subtb <- out_full_tb[out_full_tb[,1] == chr,]
  	 subtb[,1] = chr_no
  	 chr_no = chr_no +1
 	 new_tb <- rbind(new_tb, subtb)
	}
	out_file <- paste("data/indels_",S_other,".bed",sep="")
	write.table(new_tb,  out_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE )

}


## similarly extract SNPs in bed format (note that it should be 0-based for bedtools)
SNP_bed <-function(S_other)
{
	chrs <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

	## note: need to change chromosome name to numbers! (which is the how it is encoded in vcf files)

	## read the chr size file; to replace the last line for each chr
	chr_sizes <- read.table("/net/harris/vol1/data/nobackup/S288C_ref/sacCer3.chrom.sizes", stringsAsFactors=FALSE)


	out_full_tb <- NULL
	out_full_tb_no <- NULL ## chromosome in numbers
	chr_no =1
	for (chr in chrs)
	{
 	 ## read in the SNP_indel diff table per chr
 	 file_name <- paste("data/Scer_", S_other, ".sorted_",chr,"_diff_tb.txt",sep="")
  	 diff_tb<- read.table(file_name, stringsAsFactors=FALSE)
 	 SNP_out_tb <- diff_tb[diff_tb[,3]=="SNP",c(1,2,2,5)] ## add alt SNP (in Spara) into the bed file
  	
  	 out_full_tb <- rbind(out_full_tb, SNP_out_tb)
	}
	
	####### important! The bed file used in bedtools and UCSC is 0-based, and ends with 1-based
	## what my code generates for begin and end are all 1-based
	## for exaple the first base: 
	## chr1   0        1    first_base

	out_full_tb[,2] <- as.numeric(out_full_tb[,2]) -1
	
	## replace chrI ..  chrXI to 1..16 and output in bed format (to extract subset of SNPs using bedtools)
	chr_no=1
	new_tb <- NULL
	for (chr in chrs)
	{
 	 subtb <- out_full_tb[out_full_tb[,1] == chr,]
  	 subtb[,1] = chr_no
  	 chr_no = chr_no +1
 	 new_tb <- rbind(new_tb, subtb)
	}
	out_file <- paste("data/SNPs_",S_other,".bed",sep="")
	write.table(new_tb,  out_file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE )

}


## merge SNPs from different S.para, only include 
merge_SNP <- function()
{
	S_para_list = c( "YPS138", "UWOPS919171", "N44", "UFRJ50816")
	
	file_name <- paste("SNPs", "CBS432", "excl_indel.bed", sep="_")
	tb<-read.table(file_name)
	tb[,"key"] <- paste(tb[,1],tb[,3],sep="_") ## merge chr and postion (use 1-based)
	merged_tb <- tb [,c("V4","key")]
	names(merged_tb) <- c("CBS432","key")
	
	for (S_para in S_para_list)
	{
		file_name <- paste("SNPs", S_para, "excl_indel.bed", sep="_")
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
	write.table(shared_alt_allele, "merged_shared_alt_SNP_para.txt",sep="\t",quote=FALSE, row.names=F)
	
	# extract unique values (alt alleles) for each site
	alt_allele <- apply(shared_alt_allele[,c(2,3,4,5,6)], 1, function(x) unique(x[!is.na(x)]))
	library(reshape2)
	alt_allele_tb <- colsplit(shared_alt_allele[,1],pattern="_", names= c("chr","pos"))
	alt_allele_tb[,"allele"] <- alt_allele
	write.table(alt_allele_tb, "Scer_para_consol_SNPs.txt",sep="\t",quote=FALSE, row.names=F, col.names=F) ## no headers
	
	
	## for the rest of SNPs, except singletons, will output to filter out (because cannot infer ancestral state)
	
	## the table of SNPs which exclude the ones which we can call ancestral state reliably
	rest_tb <- merged_tb [!(merged_tb$key %in% shared_alt_allele$key),] 
	## get the list where there is one ALT allele called among the 5 S.para (will exclude those from filtering)
	singleton_S_para <-  merged_tb [rowSums(is.na(merged_tb)) ==4 ,]
	## the remaining SNPs are the ones that are most likely segregating in the ancestral state; will need 
	## to filter out  
	SNPs_tb_to_filter <-  rest_tb [!(rest_tb$key %in% singleton_S_para$key),]
	
	## convert the position into bed format 
	library(reshape2)
	SNPs_tb_to_filter_bed <- colsplit(SNPs_tb_to_filter[,"key"],"_",names=c("chr","pos"))
	SNPs_tb_to_filter_bed["begin"] <- SNPs_tb_to_filter_bed["pos"]-1
	write.table (SNPs_tb_to_filter_bed[,c("chr","begin","pos")], "S_para_SNP_to_filter.bed",quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}


slist <- c( "CBS432","YPS138", "UWOPS919171", "N44", "UFRJ50816")
for (strain in slist)
{
	indel_bed(strain)
}
