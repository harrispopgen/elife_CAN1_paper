## extract genes considered 

import os

from fasta_basic import *

if __name__ == '__main__':
	with open("../data/gene_in_anc.txt") as f:
			genes = [line.rstrip() for line in f]
    	
	fa_file = "/net/harris/vol1/project/yeast_1002/nobackup/S_cere_ref/archive/S288C_reference_genome_R64-1-1_20110203/orf_coding_all_R64-1-1_20110203.fasta"
	extract_genes2(fa_file, genes, "../data/orf_in_ace.fasta")

