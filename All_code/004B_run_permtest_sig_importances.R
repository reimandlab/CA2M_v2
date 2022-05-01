#Get cohort number, signature index, and block as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])
sig_index = as.integer(args[2])
block = as.integer(args[3])

source("000_HEADER.R")

#Get projects codes for PCAWG cohorts
project_codes = gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)]

#Select project code
project_code = project_codes[cohort_index]

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load in signature binned mutation counts data table for project
Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), project_code, ".csv"))[,-c(1,2)]

#Keep only sigs with > 10,000 mutations
Sigs = Sigs[,.SD, .SDcols = which(colSums(Sigs)>10000)]
                   
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

#Select signature based on input
signature = colnames(Sigs)[sig_index]
                   
if(!is.na(signature)){

    #Output for RF model is binned mutation counts for specific signature in cohort
    Output = Sigs[[signature]]

	#Remove hypermutated regions if cohort is lymphoma or leukemia
	if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
		Preds = Preds[-c(292, 2438)]
		Output = Output[-c(292, 2438)]}
	
	get_RF_dt = function(seed){
		print(seed)
		set.seed(seed)
		
		#Permute output based on random seed for permutation test
		output_scrambled = sample(as.numeric(Output))

		#Train randomForest
		rf = randomForest(x = Preds, 
						  y = output_scrambled, 
						  keep.forest = T, 
						  ntree = 1000, 
						  do.trace = F, 
						  importance = T)

		return(as.numeric(importance(rf, type = 1, scale = F)))}
	
	#Blocks for running 100 of 1000 permutations at a time
	blocks = seq(100, 1000, 100)           
	
	#Repeat permutation test 100 times
	RF_results = as.data.table(do.call("rbind.data.frame", 
									   mclapply(seq(blocks[block]-99, blocks[block]),
									   get_RF_dt, mc.cores = 10))) 
	colnames(RF_results) = colnames(Preds)
	
	#Save results
	fwrite(RF_results, paste0(pff("data/004B_Feature_importance_permutations/"),"results_",cohort_index, "_", sig_index, "_", block, ".csv"))

	print("Finished")
}else{
	file.create(paste0(pff("data/004B_Feature_importance_permutations/"),"results_",cohort_index, "_", sig_index, "_", block, ".csv"))
}
                     