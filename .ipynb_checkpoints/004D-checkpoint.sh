#!/bin/bash
for i in `seq 1 23`; do qsub -cwd -b y -P reimandlab -N shap$i -l h_vmem=5G,h_rt=3:0:0:0 -pe smp 6 -now no -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/004D_run_mutsig_RFSHAP.R $i"; done