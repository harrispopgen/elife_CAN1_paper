## first use R code to generate indels for each S.para

source loadBedtools

cat data/indels_CBS432.bed   data/indels_YPS138.bed  \
	data/indels_UWOPS919171.bed  data/indels_N44.bed  \
	data/indels_UFRJ50816.bed > data/all_indels_Spara.bed

sort -k 1,1 -k2,2n  data/all_indels_Spara.bed > data/all_indels_Spara_sorted.bed
bedtools merge -i data/all_indels_Spara_sorted.bed  > data/all_indels_Spara_merge.bed
