#!/bin/bash
for (( i = 1; i <= 17; i++ ))
do  
    for (( b = 1; b <= 10; b++ ))
    do  
        qsub -N boot-$i-$b -cwd -b y -P reimandlab  -l h_vmem=3G,h_rt=0:72:0:0 -pe smp 10 -o /dev/null -e /dev/null "source /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/etc/profile.d/conda.sh; conda activate r_env; Rscript /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/003C_RF_bootstrap_importances.R $i $b" 
    done 
done
