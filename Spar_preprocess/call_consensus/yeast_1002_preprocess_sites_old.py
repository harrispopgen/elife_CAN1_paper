## main file to preprocess
## 1. generate table from outgroup (S.paradoxus), for multiple S_para

import sys
from extract_diffs_w_chr import *

#import pdb

if __name__ == '__main__':
    print ('extracting diffs from S.cer and S.para...')
    in_dir = '../multi_para/data/'
    S_para_list=["CBS432","N44","UFRJ50816","UWOPS919171","YPS138"]
    for S_para in S_para_list:
    #pdb.set_trace()
        infile = "Scer_"+  S_para + ".sorted.net.axt"
        out_dir = 'data/'
        ## write diff per chr
        extract_diffs(in_dir, infile, out_dir)
