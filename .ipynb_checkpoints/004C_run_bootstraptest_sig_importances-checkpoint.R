#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])
sig_index = as.integer(args[2])
block = as.integer(args[3])

source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

project_codes = gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)]

project_code = project_codes[cohort_index]

print(project_code)

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load in signature data table
Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), project_code, ".csv"))[,-c(1,2)]

#Keep only sigs with > 10,000 mutations
Sigs = Sigs[,.SD, .SDcols = which(colSums(Sigs)>10000)]
                   
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}


signature = colnames(Sigs)[sig_index]

print(signature) 
                   
if(!is.na(signature)){

	#Process MAF into RMV
    Output = Sigs[[signature]]


	if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
		Preds = Preds[-c(292, 2438)]
		Output = Output[-c(292, 2438)]}
    
	get_RF_dt = function(seed){
		print(seed)
		set.seed(seed)

		#Convert input and output to data.table
		bootstrap_sample = sample(nrow(Preds), replace = T)
		input_bootstrap = Preds[bootstrap_sample]
		output_bootstrap = Output[bootstrap_sample]

		#Train randomForest
		rf = randomForest(x = input_bootstrap, 
						  y = output_bootstrap, 
						  keep.forest = T, 
						  ntree = 1000, 
						  do.trace = F, 
						  importance = T)

		return(as.numeric(importance(rf, type = 1, scale = F)))}

	blocks = seq(10, 1000, 10)  
	           
	RF_results = as.data.table(do.call("rbind.data.frame", 
									   lapply(seq(blocks[block]-9, blocks[block]),
									   get_RF_dt))) 
	colnames(RF_results) = colnames(Preds)

	fwrite(RF_results, paste0(pff("/data/004C_Sig_importance_bootstrap/"),"results_",cohort_index, "_", sig_index, "_", block, ".csv"))

	print("Finished")
}else{
	file.create(paste0(pff("/data/004C_Sig_importance_bootstrap/"),"results_",cohort_index, "_", sig_index, "_", block, ".csv"))
}
                     