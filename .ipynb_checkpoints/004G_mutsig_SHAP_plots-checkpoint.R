source("000_HEADER.R")

#Load in required python packages to load and plot SHAP values using reticulate
library(reticulate)
shap = import("shap")
plt = import('matplotlib.pyplot')
backend=import("matplotlib.backends.backend_pdf")

#Get paths to predictor importance values for signature-based analysis
importance_paths = list.files(pff("data/004A_Sig_RF_Results"), 
					   full.names = T)
imp_cohort_names = unlist(lapply(list.files(pff("data/004A_Sig_RF_Results")),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))
						 
#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           

#Load predictor matching list to PCAWG cohorts
matching_list = readRDS(pff("data/001J_Matching_predictors_dictionary.RDS"))

#Get paths of SHAP results
SHAP_paths = list.files(pff("data/004D_Sig_RF_SHAPvalues"), full.name = T)    							 
SHAP_cohort_names = list.files(pff("data/004D_Sig_RF_SHAPvalues"))							 

#Load in preds
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Get plot data for cancer type							 
plot_SHAP_ct = function(cohort_name, top_N_predictors){
	
	#Load in predictor importance values for cohort
	importances = fread(importance_paths[which(imp_cohort_names == cohort_name)])
    
	#Get predictor names
	predictors = make.unique(importances[[1]][-1])
    
	#Get signatures available
    Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), cohort_name, ".csv"))[,-c(1,2)]

	#Keep only signatures with > 20,000 mutations and ALL mutations
	Sigs_original = colnames(Sigs[,.SD, .SDcols = which(colSums(Sigs)>10000)])
	Sigs_tokeep = colnames(Sigs[,.SD, .SDcols = which(colSums(Sigs)>20000)])		
	
	#Remove hypermutated windows if cohort is lymphoma or leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
    }	
	
	#Get paths to signature-based SHAP value results
    SHAP_paths_ct = list.files(SHAP_paths[which(SHAP_cohort_names == cohort_name)], full.name = T)
    SHAP_sigs_ct = unlist(lapply(list.files(SHAP_paths[which(SHAP_cohort_names == cohort_name)]),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))
	
	#Get plot data for signature							 
	plot_SHAP_sig = function(signature){
		
		#Get importances
		data = importances[[signature]][-1]

		#Get indices for cohort and signature
		cohort_index = which(gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)] == cohort_name)
		sig_index = which(Sigs_original == signature)		
		
		#Get permutation test results										
		permutation_paths = list.files(pff("/data/004B_Feature_importance_permutations/"), 
									full.names = T)
		permutation_splits = strsplit(list.files(pff("/data/004B_Feature_importance_permutations/")), 
								   split = "_")
		permutation_files_to_keep = unlist(lapply(permutation_splits, function(x) (x[2] == cohort_index) & (x[3] == sig_index)))

		permutation_dt = as.data.table(do.call("rbind.data.frame", 
											 lapply(permutation_paths[which(permutation_files_to_keep)], fread)))
												  
		#Get bootstrap test results
		bootstrap_paths = list.files(pff("/data/004C_Sig_importance_bootstrap/"), 
									full.names = T)
		bootstrap_splits = strsplit(list.files(pff("/data/004C_Sig_importance_bootstrap/")), 
								   split = "_")
		bootstrap_files_to_keep = unlist(lapply(bootstrap_splits, function(x) (x[2] == cohort_index) & (x[3] == sig_index)))

		bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
											 lapply(bootstrap_paths[which(bootstrap_files_to_keep)], fread)))				
		#Get p-values of predictors from permutation test
		predictor_pvals = unlist(lapply(1:length(predictors), 
										function(x)  sum(data[x] < permutation_dt[[x]])/1000))
	
		#Get significance of predictors 								
		significance = ifelse(predictor_pvals == 0, "*", "")
		names(significance) = predictors
		
		#Keep top N significant predictors																
		top_N_predictors = min(top_N_predictors, sum(significance == "*"))
								
		#Get the importance mean from the boostrap test for all significant predictors								
		sig_predictor_means = unlist(lapply(which(significance == "*"), 
											function(x) mean(Bootstrap_dt[[x]])))                                   
        
		#Get the names of the top predictors									
		top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
														 decreasing = T)][1:top_N_predictors])
		
		#Get the descriptions for the top predictors
		top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_predictors, predictor_descriptions$Predictor_names)]	

		#Load in SHAP values
		SHAP_dt = fread(SHAP_paths_ct[which(SHAP_sigs_ct == signature)])

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
		plt$title(paste(cohort_name, signature))
		pdf$savefig(fig, bbox_inches = "tight")}
											
	lapply(Sigs_tokeep, plot_SHAP_sig)}

#Analyze core 17 cancer types											
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
											
#Save SHAP plots											
plt$close()                           											
pdf = backend$PdfPages(pff("data/004G_mutsig_SHAP_plots.pdf"))                          
lapply(cancer_types_to_keep, plot_SHAP_ct, 5)     
pdf$close()   											