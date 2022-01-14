source("000_HEADER.R")

#Load in required python packages to load and plot SHAP values using reticulate
library(reticulate)
shap = import("shap")
plt = import('matplotlib.pyplot')
backend=import("matplotlib.backends.backend_pdf")

#Load in top predictors
importance = fread(pff("data/003A_feature_importances_allpreds.csv"))

#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           

#Load predictor matching list to PCAWG cohorts
matching_list = readRDS(pff("data/001J_Matching_predictors_dictionary.RDS"))

#Load in permutation test paths
p_val_paths = list.files(pff("data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 
#Load in bootstrap test paths                                 
bootstrap_paths = list.files(pff("data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types = unlist(lapply(list.files(pff("data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
      
#Load in SHAP dataset paths                                     
SHAP_paths = list.files(pff("data/003F_RF_SHAP_Results/"), full.name = T)
SHAP_cancer_types = unlist(lapply(list.files(pff("data/003F_RF_SHAP_Results/")), 
								function(x) unlist(strsplit(x, split = ".csv"))[1]))

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Function to plot SHAP violin plots of importances for top N predictors for specific cohort       
plot_SHAP_plot = function(cohort_index, top_N_predictors){
   	
	#Get cohort name
	print(cohort_index)
	cohort_name = colnames(importance)[-1][cohort_index]
	
	#Get predictor names
	predictors = importance[[1]]
	
	#Get predictor importances
    data = importance[[cohort_index + 1]]
	
	#Get predictor importances from permutation test
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
	#Get permutation test p-values for each predictor (fraction of times permuted importance is more significant than actual importance)
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))
	
	#Predictor are significant if they are more important all 1000 permuted values
    significance = ifelse(predictor_pvals == 0, "*", "")
    names(significance) = predictors
									
	#Keep top N significant predictors								
	top_N_predictors = min(top_N_predictors, sum(significance == "*"))
	
	#Get predictor importances from bootstrap test
    Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
														   full.name = T), fread)))
	#Get the importance mean from the boostrap test for all significant predictors								
    sig_predictor_means = unlist(lapply(which(significance == "*"), 
										function(x) mean(Bootstrap_dt[[x]])))                                   
    
	#Rank significant predictors by bootstrap mean
	top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
													 decreasing = T)][1:top_N_predictors])
										
	#Get predictor importance means and sd from bootstrap experiment for each significant predictor                     
    top_predictor_means = as.numeric(sig_predictor_means[order(sig_predictor_means, 
															   decreasing = T)][1:top_N_predictors])
    top_predictor_SDs = unlist(lapply(top_predictors, function(x) sd(Bootstrap_dt[[x]])))								

    #Get predictor categories from supplementary file
    top_predictor_categories = paste0(predictor_descriptions$Predictor_categories[match(top_predictors, 
																	   predictor_descriptions$Predictor_names)], " CA")
	top_predictor_categories = ifelse(top_predictor_categories=="RT CA", 
									  "Replication timing", top_predictor_categories)
	
	#Check if each predictor is matching the tissue type of cohort
    top_predictor_matching = ifelse(top_predictors %in% matching_list[[cohort_name]], 
									"Matching", 
									"Non-Matching")
	
	#Get colour category of each predictor based on predictor category and matching status								  
    top_predictor_fill = unlist(lapply(1:length(top_predictors), 
									   function(x) ifelse(top_predictor_matching[x]=="Matching", 
														  paste(top_predictor_categories[x], "from matching tissue", sep = " "), 
														  paste(top_predictor_categories[x], " from non-matching tissue", sep = ""))))	
	#Get predictor descriptions								   
	top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_predictors, predictor_descriptions$Predictor_names)]	
	
    #Load in SHAP values
    SHAP_dt = fread(SHAP_paths[which(SHAP_cancer_types == cohort_name)])
										
    #Keep SHAP values for top predictors
    SHAP_dt_top = SHAP_dt[,.SD,.SDcols = top_predictors]
    colnames(SHAP_dt_top) = top_predictor_descriptions
    Input_top = Preds[,.SD,.SDcols = top_predictors]                                  
    colnames(Input_top) = top_predictor_descriptions
    
    #Plot SHAP summary plot                                  
    fig = plt$figure(figsize = c(3, 3))
    shap$summary_plot(data.matrix(SHAP_dt_top), 
					  Input_top, 
					  show = F,
					  sort = F, 
					  max_display = 5L)
    plt$title(cohort_name)
    pdf$savefig(fig, bbox_inches = "tight")                                  
                                      
}
                                      
#Get SHAP plots for all cancer types                                      
plt$close()                           
pdf = backend$PdfPages(pff("data/003G_toppred_SHAP_plots.pdf"))                          
lapply(seq(26)[-8], plot_SHAP_plot, 5)     
pdf$close()                                                               