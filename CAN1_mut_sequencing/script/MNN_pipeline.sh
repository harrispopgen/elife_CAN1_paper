## pipeline for calling MNMs

##find candidate MNMs, with pairs of variants that are close to each other (within 10bp), and have similar allele frequency in each pool.

Rscript get_MNM_candidate.R

## validate candidate MNMs by looking at the aligned reads in the sam/bam file, and calculate the percentage 
## of occurrence of the two mutations (tried with 2 different MQ cutoff)

bash validate_MNM_samples.sh

## consolidate MNMs by merging multiple MNMs that have shared mutations in the same pool 

Rscript consolidate_MNM.R



