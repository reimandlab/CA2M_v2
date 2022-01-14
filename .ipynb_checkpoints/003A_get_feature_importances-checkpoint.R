source("000_HEADER.R")

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

#Get predictor importances and model accuracy for random forest model
get_importances = function(cohort_index){
	
	#Get cohort name using cohort index
    project_code = colnames(mutation_counts_dt)[cohort_index] 
	
	#Output of RF model is the binned mutation counts of cohort of interest
    output = mutation_counts_dt[[cohort_index]]

    #Remove hypermutated windows if cohort is lymphoma or leukemia
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
        output = output[-c(292, 2438)]
    }
    
    #Train randomForest model
    rf = randomForest(x = Preds, 
					  y = output, 
					  keep.forest = T,
					  ntree = 1000, 
					  do.trace = F,
					  importance = T)
    
	#Get predictor importances  (incMSE)
    importances = as.numeric(importance(rf, type = 1, scale = F))
    
	#Get model accuracy as R2 score
    R2 = (cor(rf$predicted, output))**2
    
	#Convert to adjusted R2
    adj_R2 = adjust_R2(R2, nrow(Preds), ncol(Preds))
    
	#Combine into list
    results = c(adj_R2, importances)
    
    return(results)
}

#Get importances and accuracies for all cohorts
result_dt = as.data.table(do.call("cbind.data.frame", 
								  mclapply(seq(26), get_importances, mc.cores = 4)))
colnames(result_dt) = colnames(mutation_counts_dt)
                   
result_dt = as.data.table(cbind.data.frame(Pred = c("Adj_R2", colnames(Preds)), 
										   result_dt))                   

#Save results
fwrite(result_dt, pff("/data/003A_feature_importances_allpreds.csv"))

#Get top 5 features
# apply(RF_importances[,-1], 2, function(x) RF_importances[,1][order(x, decreasing = T)[1:5]] )


                                                                 
