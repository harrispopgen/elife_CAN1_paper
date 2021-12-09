## based on Kelley's code of "extract_diffs.py"
## added chromosome field; modified for deletion in alignment on target or source strand
## generated one file; with Indel meaning the SNPs that needs to be filtered out!
## for indels, rather than output the length, output the begin and end positions

##axt file: coordinates begin with 1, ranges are both inclusive

import pdb

## note: the infile ends with .net.axt
## this axt file contains all the chromosome alignment info
def extract_diffs (input_dir, input_file_name, out_dir):

    input_file = input_dir + input_file_name
    infile = open(input_file,"r")
    lines=infile.readlines()
    infile.close()

    line_ind=0
    while not lines[line_ind][0]=='0':
        line_ind+=1

    ## add chromosome
    output=""
    ## do not output header
    #output='Chr Pos SNP/Indel Target Query\n' ## this is the output with table SNP with ancestral and derived allele
    ## for indel : chr; begin; indel; begin; end; (begin and end are both inclusive; repeat begin twice, to have the same number of fields as SNP)

    #pdb.set_trace()
    last_bin_start=0
    last_end=0
    last_chr = lines[line_ind].split(' ')[1]

    while line_ind<len(lines):
        s=lines[line_ind].split(' ')
        ## current chromosome
        chr = s[1]

        if last_chr!= chr:
            ## write to file
            output+= last_chr + " " + str(last_end+1)+' Indel ' + str(last_end+1) + ' 100000000\n'

            ## exclude axt
            outfile_name =  out_dir + input_file_name[:-8]+ "_" + last_chr + "_diff_tb.txt"

            outfile=open(outfile_name,'w')
            outfile.write(output)
            outfile.close()

            output ="" ##clear out
            last_chr = chr

    #    print s
        start, end = int(s[2]), int(s[3])
        if start - last_end > 1: ## if there is a gap between two alignment, add to filter
            output+= chr + " " + str(last_end+1)+' Indel '  + str(last_end+1) + " "   +str(start-1)+'\n'
        last_end=end
        if start>last_bin_start+10**6:
            # print start
            last_bin_start=start
        target=lines[line_ind+1].upper().strip('\n')
        query=lines[line_ind+2].upper().strip('\n')
        char_ind=0
        while char_ind<len(target):
            found_indel=False
            t_allele=''
            q_allele=''
    #        print target[char_ind], query[char_ind]
            if target[char_ind]==query[char_ind] or target[char_ind]=='N' or query[char_ind]=='N':
                char_ind+=1
            elif target[char_ind] in 'ACGT' and query[char_ind] in 'ACGT':
                found_match=False
                t_allele=target[char_ind]
                q_allele=query[char_ind]
                char_ind+=1
                while not found_match and not found_indel and char_ind<len(target):
                    if target[char_ind]==query[char_ind]:
                        found_match=True
                    elif target[char_ind]=='-' or query[char_ind]=='-':
                        found_indel=True
                    else:
                        t_allele+=target[char_ind]
                        q_allele+=query[char_ind]
                        char_ind+=1
                if not found_indel:
                    for i in range(len(t_allele)):
                        j=char_ind-len(t_allele)+i
                        output+= chr + " "  + str(start+j)+' SNP '+target[j]+' '+query[j]+'\n'
                else : ## if there is indel right after SNP, mark the whole region to the filtered table
                    output += chr + " "  + str(start + char_ind-len(t_allele)) + ' Indel ' + str(start + char_ind-len(t_allele)) + " "+ str(start + char_ind-1)+'\n'

            if char_ind<len(target) and (target[char_ind]=='-' or query[char_ind]=='-'):
                indel_start = start+char_ind
                len_indel=0
                if (target[char_ind]=='-'): ## missing base is on target
                    flag=1
                else:
                    flag=2
                len_SNPs_follow=0
                while char_ind<len(target) and not target[char_ind]==query[char_ind]:
                    if (target[char_ind]=='-' or query[char_ind]=='-'):
                        end = "del"
                    else:
                        len_SNPs_follow +=1
                    if target[char_ind]=='-':
                        start-=1
                    len_indel+=1
                    char_ind+=1
                if (flag==1 and len_SNPs_follow!=0):  ## if missing base is on target, then mask SNPs to follow; if missing is on query, mask both del and SNPs
                    output+= chr + " "  + str(indel_start) + ' Indel ' + str(indel_start) + " " + str(indel_start + len_SNPs_follow-1)+'\n'
                elif flag ==2:
                    output+= chr + " "  + str(indel_start) + ' Indel ' + str(indel_start) + " " + str(indel_start + len_indel-1)+'\n'

        line_ind+=4

    ## write the last chr to file
    output+= chr + " " + str(last_end+1)+' Indel ' + str(last_end+1) + ' 100000000\n'

    ## exclude axt
    outfile_name =  out_dir + input_file_name[:-8]+ "_" + chr + "_diff_tb.txt"

    outfile=open(outfile_name,'w')
    outfile.write(output)
    outfile.close()


# if __name__ == '__main__':
#     import sys
#     chrom = int(sys.argv[1])
#     print 'extracted diffs for chromosome %i' % chrom
#     infile = '../data/hg19_chimp_align/chr'+ str(chrom) + '.hg19.panTro4.net.axt'
#     extract_diffs(infile)
