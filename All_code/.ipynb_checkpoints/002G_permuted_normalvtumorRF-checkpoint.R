#Get cohort number as argument to loop through
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("000_HEADER.R")

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load in all predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Core cancer types to keep
cancer_types_to_keep = c("PANCAN", "Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
#Load megebase-scale binned mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 

#Get cohort name
project_code = cancer_types_to_keep[cohort_index]

#Response of model will be mutation counts from cohort of interest
output = mutation_counts_dt[[project_code]]                                    

#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           
tumor_predictor_descriptions = predictor_descriptions[`Predictor_categories` == "Primary cancer"]
normal_predictor_descriptions = predictor_descriptions[`Predictor_categories` == "Normal tissue"]

#Get total number of priamry cancer predictors matching to core cancer types
tumor_counts = sapply(cancer_types_to_keep, function(x)
	length(grep(x, tumor_predictor_descriptions$Predictor_matching_PCAWG)))

#Get total number of normal cell predictors matching to core cancer types
normal_counts = sapply(cancer_types_to_keep, function(x)
	length(grep(x, normal_predictor_descriptions$Predictor_matching_PCAWG)))
					  
#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

#Function to permute predictors based on minimum number of predictors in either tumor or normal set
get_predictors = function(cancer_type, seed){
	
	set.seed(seed)
	
	#Get minimum number of predictors for cancer type
	min_preds = min(tumor_counts[cancer_type], normal_counts[cancer_type])
	
	#Permute minimum number of both tumor and normal predictors
	tumor_preds_sampled = sample(tumor_predictor_descriptions$Predictor_names[grep(cancer_type,tumor_predictor_descriptions$Predictor_matching_PCAWG)], min_preds)
	normal_preds_sampled = sample(normal_predictor_descriptions$Predictor_names[grep(cancer_type, normal_predictor_descriptions$Predictor_matching_PCAWG)], min_preds)
	
	#Return predictor names
	return(list(tumor_preds_sampled, normal_preds_sampled))}
	
	
	
#Run random forest on tumor and normal predictors					   
run_RF = function(seed, type){
	
	###Type: 1 for tumor preds and 2 for normal preds
	
	#Get predictor names sampled and matched for each core cancer_type using seed
	pred_names_to_keep = unique(unlist(lapply(cancer_types_to_keep, function(x) 
		get_predictors(x, seed)[[type]])))

	#Get replication timing predictor names
	RT_pred_names = predictor_descriptions[`Predictor_categories` == "RT"]$Predictor_names
											  
	Preds_new = Preds[ , .SD, .SDcols = c(pred_names_to_keep , RT_pred_names)]											  
	
	#Remove hypermutated windows if cohort is lymphoma or leukemia
    if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds_new = Preds_new[-c(292, 2438)]
        output = output[-c(292, 2438)]
    }
    
    #Train randomForest model
    rf = randomForest(x = Preds_new, 
					  y = output, 
					  keep.forest = T,
					  ntree = 1000, 
					  do.trace = F,
					  importance = F)
    
	#Get model accuracy as R2 score
    R2 = (cor(rf$predicted, output))**2
    
	#Convert to adjusted R2
    adj_R2 = adjust_R2(R2, nrow(Preds_new), ncol(Preds_new))
    
    return(adj_R2)
}

#Get importances and accuracies for tumor trained models
tumor_result_dt = unlist(mclapply(1:1000, run_RF, 1, mc.cores = 8))
              
#Save tumor trained model results
saveRDS(tumor_result_dt, pff(paste0("/data/002G_permuted_normalvtumor_results/Tumor_results/", project_code, ".RDS")))
											 
#Get importances and accuracies for normal trained models
normal_result_dt = unlist(mclapply(1:1000, run_RF, 2, mc.cores = 8))
              
#Save tumor trained model results
saveRDS(normal_result_dt, pff(paste0("/data/002G_permuted_normalvtumor_results/Normal_results/", project_code, ".RDS")))
			