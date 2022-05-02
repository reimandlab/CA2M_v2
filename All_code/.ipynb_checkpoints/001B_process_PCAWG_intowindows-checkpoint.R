source("000_HEADER.R")

source(pff("/bin/999_process_data.R"))

#Load in hg38 PCAWG MAF file
PCAWG_mutations_dt_hg38 = fread(pff("data/001A_PCAWG_mutations_hg38.csv"))

#Define cohort sizes for PCAWG
PCAWG_sample_table = PCAWG_mutations_dt_hg38[,.(Cohort_size = uniqueN(Donor_ID)), by = Project_Code]

#Keep cohorts with >30 patients
cohorts_to_keep = sort(PCAWG_sample_table[Cohort_size>30]$Project_Code)

#Convert SNV's to genomic window-based mutation counts
process_muts = function(window_size){

    get_mutation_counts = function(cohort){

        mut = PCAWG_mutations_dt_hg38[Project_Code == cohort] #Keep mutations for specific cohort

        mut_counts = create_RMV(mut, window_size)[[3]] #Convert to mutation counts based on window size

    }
	
	#Get cohort-specific window counts
    mutation_counts_dt = as.data.table(do.call("cbind",lapply(cohorts_to_keep, get_mutation_counts)))
	
	#Add pan-cancer window counts
    mutation_counts_dt = as.data.table(cbind.data.frame(create_RMV(PCAWG_mutations_dt_hg38, window_size), mutation_counts_dt))
    colnames(mutation_counts_dt) = c("chr", "start", "PANCAN", cohorts_to_keep)

    #Remove windows with less than or equal to 80% mappability
    mutation_counts_dt_mappable = filter_windows_with_UMAP(mutation_counts_dt, window_size)
    
    return(mutation_counts_dt_mappable)}

#Save mutation counts files for different window sizes
fwrite(process_muts(1000000), pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))
fwrite(process_muts(100000), pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))