#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("000_HEADER.R")

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
                
#Create function to run RF
run_RF = function(signature){

	print(signature)

    #Output for RF model is binned mutation counts for specific signature in cohort
    Output = Sigs[[signature]]
    
	#Remove hypermutated regions if cohort is lymphoma or leukemia
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
        Output = Output[-c(292, 2438)]}
    
    #Train randomForest
    rf = randomForest(x = Preds, 
					  y = Output, 
					  keep.forest = T, 
					  ntree = 1000, 
					  do.trace = F, 
					  importance = T)
    
	#Get predictor importances as IncMSE
    importances = as.numeric(importance(rf, type = 1, scale = F))
    
	#Get R2 accuracy
    R2 = (cor(rf$predicted, Output))**2
    
	#Adjust R2
    adj_R2 = adjust_R2(R2, nrow(Preds), ncol(Preds))
    
    results = c(adj_R2, importances)
    
    return(results)}

#Run for all signatures
result_dt = as.data.table(do.call("cbind.data.frame", 
								  mclapply(colnames(Sigs), run_RF, mc.cores = 4)))
colnames(result_dt) = colnames(Sigs)
                   
result_dt = as.data.table(cbind.data.frame(Pred = c("Adj_R2", colnames(Preds)), 
										   result_dt))                   
                   
print("Finished")
  
#Save results
fwrite(result_dt, 
	   paste0(pff("data/004A_Sig_RF_Results/"), project_code, ".csv"))