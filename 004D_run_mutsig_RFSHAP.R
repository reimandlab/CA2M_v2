#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
source(pff("/bin/999_run_RF_withSHAP.R"))

project_codes = gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)]

project_code = project_codes[cohort_index]

print(project_code)

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load in signature data table
Sigs_dt = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), project_code, ".csv"))[,-c(1,2)]

#Keep only sigs with > 10,000 mutations
Sigs = colnames(Sigs_dt[,.SD, .SDcols = which(colSums(Sigs_dt)>10000)])
                   
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}
                  
dir.create(paste0(pff("data/004D_Sig_RF_SHAPvalues/"), project_code, "/"))
                   
run_sig_RF = function(signature){
    
    #Process MAF into RMV
    Output = Sigs_dt[[signature]]
    
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292,2438)]
        Output = Output[-c(292,2438)]}
    
    SHAP_values = run_RF_and_get_SHAP(Preds, Output)

    fwrite(SHAP_values, paste0(pff("data/004D_Sig_RF_SHAPvalues/"), 
							   project_code, "/", signature, ".csv"))

}

mclapply(Sigs, run_sig_RF, mc.cores = 6)

print("Finished")}

					 