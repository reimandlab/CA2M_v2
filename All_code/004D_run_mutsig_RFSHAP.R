#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("000_HEADER.R")

source(pff("/bin/999_run_RF_withSHAP.R"))

#Get projects codes for PCAWG cohorts
project_codes = gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)]

#Select project code
project_code = project_codes[cohort_index]

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load in signature binned mutation counts data table for project
Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), project_code, ".csv"))[,-c(1,2)]

#Keep only signatures with > 10,000 mutations in cohort
Sigs = Sigs[,.SD, .SDcols = which(colSums(Sigs)>10000)]
                   
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

#Create folder for results
dir.create(paste0(pff("data/004D_Sig_RF_SHAPvalues/"), project_code, "/"))

#Function to get SHAP values for predictors when predicting signature-based mutation counts
run_sig_RF = function(signature){
    
    #Process MAF into RMV
    Output = Sigs_dt[[signature]]
    
    #Remove hypermutated regions if cohort is lymphoma or leukemia
	if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292,2438)]
        Output = Output[-c(292,2438)]}
    
	#Get SHAP values
    SHAP_values = run_RF_and_get_SHAP(Preds, Output)
	
	#Save SHAP values
    fwrite(SHAP_values, paste0(pff("data/004D_Sig_RF_SHAPvalues/"), 
							   project_code, "/", signature, ".csv"))

}

#Get SHAP values for all signatures
mclapply(Sigs, run_sig_RF, mc.cores = 6)

print("Finished")}

					 