source ("chr_cvt_basic.R")

indel_bed <- read.table("data/all_indels_Spara_merge.bed", stringsAsFactors=FALSE)

indel_bed_converted <- chr_num2rom(indel_bed)

write.table(indel_bed_converted, "data/all_indels_Spara_Merge.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
