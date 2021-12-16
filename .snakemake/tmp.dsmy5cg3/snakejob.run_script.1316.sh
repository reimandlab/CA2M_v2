#!/bin/sh
# properties = {"type": "single", "rule": "run_script", "local": false, "input": ["/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/004C_run_bootstraptest_sig_importances.R"], "output": ["/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/data/004C_Sig_importance_bootstrap/results_11_12_6.csv"], "wildcards": {"variable_i": "11", "variable_s": "12", "variable_b": "6"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp", "threads": "1", "runtime": "4:0:0:0", "individual_core_memory": "15G"}, "jobid": 1316, "cluster": {}}
 cd /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin && \
/.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/envs/r_env/bin/python3.7 \
-m snakemake /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/data/004C_Sig_importance_bootstrap/results_11_12_6.csv --snakefile /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/Snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/.snakemake/tmp.dsmy5cg3 /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/004C_run_bootstraptest_sig_importances.R --latency-wait 45 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules run_script --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/envs/r_env/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/.snakemake/tmp.dsmy5cg3/1316.jobfinished || (touch /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/.snakemake/tmp.dsmy5cg3/1316.jobfailed; exit 1)

