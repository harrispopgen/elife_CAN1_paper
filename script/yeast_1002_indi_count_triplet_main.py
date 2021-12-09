## this is the final code to generate mutation spectrum for single base pair mutations
## using the same strains as used in AC=4 cutoff

from count_mutation_types_indi_new import *

if __name__ == '__main__':
    chromsomes = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

    ## new vcf
    vcf_prefix = "../data/vcf/1011_SNP_bi_mm0.8_GT_u_rm_final_no_close.chr"

    outgroup_file = "../Spar_preprocess/call_consensus/data/Scer_para_consol_SNPs.txt"

    outgroup_df = pd.read_table(outgroup_file, header = None,  sep="\t")

    ref_prefix ="/net/harris/vol1/data/Pengyao/nobackup/S288C_ref/chromosomes/chr"

    print ("chromosome 1")
    vcf_file = vcf_prefix + str(1) + ".vcf.gz"

    ref_chr = ref_prefix + chromsomes[0] + ".fa.gz"

    outgroup_tb = outgroup_df[outgroup_df[0] == ("chr" + chromsomes[0] )]

    mutaton_counts1, mutation_counts2, haplo_ids = count_mutations_triplet_chr( ref_chr, vcf_file, outgroup_tb ,1, [],[])

    for i in range(1,len(chromsomes)):
        print ("chromosome " + str (i+1))
        vcf_file = vcf_prefix + str(i+1) + ".vcf.gz"

        ref_chr = ref_prefix + chromsomes[i] + ".fa.gz"

        outgroup_tb = outgroup_df[outgroup_df[0] == ("chr" + chromsomes[i] )]

        # vcf_file = "../data/vcfs/test.vcf"
        mutaton_counts1, mutation_counts2, haplo_ids = count_mutations_triplet_chr(ref_chr, vcf_file, outgroup_tb, 0, mutaton_counts1, mutation_counts2)

    write_output2(mutaton_counts1, mutation_counts2, mutation_types, mutation_dict, haplo_ids, "../data/Scer_triplet_no_close_indi_mut_count_strain.txt")
