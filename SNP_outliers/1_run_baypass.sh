#!/bin/bash

### Process SNP data using baypass

## use poolfstat to extract the SNPs, write to a baypass file
# pooldata2genobaypass(Evol_700, writing.dir = getwd(), subsamplesize = -1)

# run this one first to get omega matrix etc
baypass_2.3/sources/i_baypass -npop 4 -gfile Evol700.genobaypass -poolsizefile Evol700.poolsize -d0yij 8 -outprefix Evol700 -nthreads 64

# then run this to look at population differences
# groups file: 0 1 0 1
baypass_2.3/sources/i_baypass  -nthreads 64 -npop 4 -gfile Evol700.genobaypass -efile groups -poolsizefile Evol700.poolsize -d0yij 8 -auxmodel -omegafile Evol700_mat_omega.out -outprefix Evol700.grps.aux
