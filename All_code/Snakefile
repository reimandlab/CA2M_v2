# OLIVER STARTS HERE

# r paths and script
r_path = '/.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/envs/r_env/bin/Rscript'

# remember that python includes 0 and goes from [start : (stop - 1)]
i = list(map(str, list(range(24))))[1:]

s = list(map(str, list(range(13))))[1:]

b = list(map(str, list(range(101))))[1:]


# create the rule that requires the output files 
rule all:
    input:
        expand("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/data/004C_Sig_importance_bootstrap/results_{variable_i}_{variable_s}_{variable_b}.csv", variable_i = i, variable_s = s, variable_b = b)


# create the qsub command
rule run_script:
    input:
        # r_path = '/.mounts/labs/reimandlab/private/users/oocsenas/Anaconda/conda/envs/r_env/bin/Rscript',
        r_script = '/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/004C_run_bootstraptest_sig_importances.R'
    output:
            '/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/data/004C_Sig_importance_bootstrap/results_{variable_i}_{variable_s}_{variable_b}.csv'
    resources:
        threads = '1',
        runtime = '4:0:0:0',
        individual_core_memory = '15G'
    wildcard_constraints:
        variable_i="[0-9]+",
        variable_s="[0-9]+",
        variable_b="[0-9]+"
    shell:
        "{r_path} {input.r_script} {wildcards.variable_i} {wildcards.variable_s} {wildcards.variable_b}"



# nohup snakemake --cluster "qsub -cwd -b y -P reimandlab -l h_vmem={resources.individual_core_memory},h_rt={resources.runtime} -pe smp {resources.threads} -o /dev/null -e /.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/snake_errors/" -j 400 -w 45 &


