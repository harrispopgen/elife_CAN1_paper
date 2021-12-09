## simulate mutations by the AC count
## input: (1) total number of simulations. (2-4) #of mutations for AC= 2, 3, 4

import random
import pdb
import pandas as pd
import numpy as np
import sys

## AC_count_vec has 3 elements, of mut counts for AC= 2, 3, 4
## only record the proportion of C->A mutations
def simulate_by_AC_CtoAratio(full_mut_tb, AC_count_vec, nsimu):
    AC_count = [2,3,4]
    nCtoA_vec = np.zeros(nsimu)

    for i in range(len(AC_count_vec)):
        #pdb.set_trace()
        AC_num = AC_count[i]
        mut_count = AC_count_vec[i]
        ## get the sub tb
        mut_tb_AC = full_mut_tb[full_mut_tb["AC"] == AC_num]
        for j in range(nsimu):
            rd_idx = random.choices(mut_tb_AC.index, k=mut_count)
            ## get the random list of mutations
            rd_mut = mut_tb_AC.loc[rd_idx]
            ## count the number of CtoA and non-CtoA mutations
            nCtoA= len(rd_mut[rd_mut["mut"] == "C_A"])
            nCtoA_vec [j] += nCtoA

    ## calculate the ratio in the end
    CtoAratio = nCtoA_vec/sum(AC_count_vec)
    return CtoAratio

if __name__ == '__main__':
    nsim = int(sys.argv[1])
    strain = sys.argv[2]
    AC2 = int(sys.argv[3])
    AC3 = int(sys.argv[4])
    AC4 = int(sys.argv[5])


    #pdb.set_trace()
    full_mut_tb = pd.read_table("../data/mut_info_AC_4_rev_comp.csv",  sep=",")
    CtoAratio = simulate_by_AC_CtoAratio(full_mut_tb, [AC2,AC3,AC4], nsim)

    outfile = "../data/simu/"  + strain + "_CA_ratio.txt"
    df = pd.DataFrame({strain: CtoAratio})
    df.to_csv(outfile,  index=False)
