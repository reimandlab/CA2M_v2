source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

source(pff("/bin/999_process_data.R"))

#Load in hg38 PCAWG file
PCAWG_mutations_dt_hg38 = fread(pff("data/PCAWG_mutations_hg38.csv"))

#Load in Cohort sizes for PCAWG
PCAWG_sample_table = PCAWG_mutations_dt_hg38[,.(Cohort_size = uniqueN(Donor_ID)), by = Project_Code]
cohorts_to_keep = sort(PCAWG_sample_table[Cohort_size>30]$Project_Code)

process_muts = function(window_size){

    get_mutation_counts = function(cohort){

        mut = PCAWG_mutations_dt_hg38[Project_Code == cohort]

        mut_counts = create_RMV(mut, window_size)[[3]]

    }

    mutation_counts_dt = as.data.table(do.call("cbind",lapply(cohorts_to_keep, get_mutation_counts)))

    mutation_counts_dt = as.data.table(cbind.data.frame(create_RMV(PCAWG_mutations_dt_hg38, window_size), mutation_counts_dt))
    colnames(mutation_counts_dt) = c("chr", "start", "PANCAN", cohorts_to_keep)

    #Remove windows with less than or equal to 80% mappability
    mutation_counts_dt_mappable = filter_windows_with_UMAP(mutation_counts_dt, window_size)
    
    return(mutation_counts_dt_mappable)}

#Save mutation counts file
fwrite(process_muts(1000000), pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))
fwrite(process_muts(100000), pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))