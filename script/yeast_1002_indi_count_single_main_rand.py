## modified on 9.10.18, which require to have outgroup file in count_mutation functions
## test whether it reproduce the same output as Kelley's code for human 1000 genomes mutation spectrum
from count_mutation_types_indi_new import *

import pdb

if __name__ == '__main__':

    chromsomes = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]

    ## new vcf
    vcf_prefix = "../data/vcf/1011_SNP_bi_mm0.8_GT_u_rm_final_no_close.chr"

    outgroup_file = "../Spar_preprocess/call_consensus/data/Scer_para_consol_SNPs.txt"

    outgroup_df = pd.read_table(outgroup_file, header = None,  sep="\t")


    print ("chromosome 1")
    vcf_file = vcf_prefix + str(1) + ".vcf.gz"

    outgroup_tb = outgroup_df[outgroup_df[0] == ("chr" + chromsomes[0] )]

    #pdb.set_trace()
    mutaton_counts1, mutation_counts2, haplo_ids = count_mutations_single_chr(vcf_file, outgroup_tb ,1, [],[], 2,0)

    for i in range(1,len(chromsomes)):
        print ("chromosome " + str (i+1))
        vcf_file = vcf_prefix + str(i+1) + ".vcf.gz"

        outgroup_tb = outgroup_df[outgroup_df[0] == ("chr" + chromsomes[i] )]

        # vcf_file = "../data/vcfs/test.vcf"
        mutaton_counts1, mutation_counts2, haplo_ids = count_mutations_single_chr(vcf_file, outgroup_tb, 0, mutaton_counts1, mutation_counts2 ,2,0)

    write_output2(mutaton_counts1, mutation_counts2, mutation_types_single, mutation_dict_single, haplo_ids, "../data/Scer_single_no_close_indi_mut_count_rand.txt" )
