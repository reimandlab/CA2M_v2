#!/bin/bash
for i in `seq 1 17`; do qsub -cwd -b y -P reimandlab -N mutsig$i -l h_vmem=20G,h_rt=0:48:0:0 -pe smp 4 -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/005P_getmutsig_100kB.R $i"; done