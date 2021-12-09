## extract all the sites that are aligned! (exclude N or indel)
## will output in bed format (so that I could use bedtools to intersect it easily)
## will use all merged indels (generated in preprocess2/) to filter out sites in the end!


##axt file: coordinates begin with 1, ranges are both inclusive

import pdb

## note: the infile ends with .net.axt
## this axt file contains all the chromosome alignment info
def extract_aln_sites (input_dir, input_file_name, out_dir):

    input_file = input_dir + input_file_name
    infile = open(input_file,"r")
    lines=infile.readlines()
    infile.close()

    line_ind=0
    while not lines[line_ind][0]=='0':
        line_ind+=1

    ## add chromosome
    #output=""
    ## do not output header

    #pdb.set_trace()

    ## exclude axt
    outfile_name =  out_dir + input_file_name[:-8]+ "_sites.bed"

    outfile=open(outfile_name,'w')

    while line_ind<len(lines):
        s=lines[line_ind].split(' ')
        ## current chromosome
        chr = s[1]

    #    print s
        start, end = int(s[2]), int(s[3])

        target=lines[line_ind+1].upper().strip('\n')
        query=lines[line_ind+2].upper().strip('\n')
        char_ind=0
        ref_gap_len=0
        output_str=""
        while char_ind<len(target):

    #        print target[char_ind], query[char_ind]
            if target[char_ind]=='N' or query[char_ind]=='N' or query[char_ind]=='-':
                char_ind+=1
            elif target[char_ind]=='-': ## a gap in the ref
                ref_gap_len +=1
                char_ind+=1
            elif target[char_ind] in 'ACGT' and query[char_ind] in 'ACGT': # output the target seq
                output_str += chr + "\t" + str(start + char_ind- ref_gap_len -1 ) + "\t" + str(start + char_ind- ref_gap_len ) + "\t" + query[char_ind] + "\n"
                char_ind+=1
            else: ## if not in any of these cases
                print ("exception!!")

        line_ind+=4
        outfile.write(output_str)

    outfile.close()
