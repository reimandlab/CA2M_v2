#Get cohort number as argument to loop through
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("000_HEADER.R")

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load megebase-scale binned mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 

#Response of model will be mutation counts from cohort of interest
Output = mutation_counts_dt[[cohort_index]]                                    
 
#Get cohort name
project_code = colnames(mutation_counts_dt)[cohort_index]

#Load in normal tissue ATAC-Seq and replication timing predictors
Normal_preds = fread(pff("data/001G_Normal_preds_1MB.csv"))[,-c(1,2)]
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]
Preds = as.data.table(cbind.data.frame(Normal_preds, RT))

#Remove hypermutated windows if cohort is lymphoma or leukemia
if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
    Preds = Preds[-c(292, 2438)]
    Output = Output[-c(292, 2438)]
}

#Run random forest with monte carlo cross-validation
RF_result = run_RF_and_get_pred_importances(Preds, Output, 1000, n_tree = 1000, cores = 16, train_split = 0.8)

#Save results
fwrite(RF_result, paste0(pff("/data/002B_normalcell_RF_Results/"), project_code, ".csv"))

print("Finished")