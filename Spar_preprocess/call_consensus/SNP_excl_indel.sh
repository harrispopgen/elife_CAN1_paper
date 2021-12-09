## exclude indels 

source loadBedtools

bedtools intersect -v -a data/Scer_CBS432.sorted_sites.bed  -b  data/all_indels_Spara_Merge.bed > data/CBS432_sites_excl_indel.bed
bedtools intersect -v -a data/Scer_YPS138.sorted_sites.bed  -b  data/all_indels_Spara_Merge.bed  > data/YPS138_sites_excl_indel.bed
bedtools intersect -v -a data/Scer_UWOPS919171.sorted_sites.bed -b  data/all_indels_Spara_Merge.bed > data/UWOPS919171_sites_excl_indel.bed
bedtools intersect -v -a data/Scer_N44.sorted_sites.bed  -b  data/all_indels_Spara_Merge.bed >  data/N44_sites_excl_indel.bed
bedtools intersect -v -a data/Scer_UFRJ50816.sorted_sites.bed  -b data/all_indels_Spara_Merge.bed > data/UFRJ50816_sites_excl_indel.bed
