#!/bin/bash
for i in `seq 1 17`; do qsub -cwd -b y -P reimandlab -N sig$i -l h_rt=500000 -l h_vmem=35g -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/005S_run_mutsig_RF_100KB.R $i"; done