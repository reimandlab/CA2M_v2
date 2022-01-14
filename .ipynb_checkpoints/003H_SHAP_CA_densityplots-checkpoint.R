source("000_HEADER.R")

#Load in required python packages to load and plot SHAP values using reticulate
library(reticulate)
shap = import("shap")
plt = import('matplotlib.pyplot')
backend=import("matplotlib.backends.backend_pdf")

#Load in top preds
importance = fread(pff("data/003A_feature_importances_allpreds.csv"))

#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           

#Load predictor matching list
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

#Get data for SHAP vs. CA plot								  
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
    Input_top = as.data.table(apply(Preds[,.SD, .SDcols = top_predictors], 2, function(x) scale(log(1+x)) ))
    colnames(Input_top) = top_predictor_descriptions  
										
    #Create tall format                                  
    SHAP_dt.m = melt(SHAP_dt_top)
    Input_dt.m = melt(Input_top)
                                      
    #Combine to final dt
    plot_dt.m = as.data.table(cbind.data.frame(variable = SHAP_dt.m$variable, 
											   Input = Input_dt.m$value, 
											   SHAP = SHAP_dt.m$value, 
                                             cohort_name = rep(cohort_name, 
																nrow(Input_dt.m))))
									
	plot_dt.m$type = predictor_descriptions$Predictor_categories[match(plot_dt.m$variable, predictor_descriptions$Predictor_descriptions)]									
    
    return(plot_dt.m)}
                                      
#Get SHAP vs. CA data for all cancer types
plot_dt = as.data.table(do.call("rbind.data.frame",lapply(seq(26)[-8], plot_SHAP_plot, 5)))

#Keep data for only core cancer types									
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")					  

plot_dt_select = plot_dt[cohort_name %in% cancer_types_to_keep]
									
#Separate data based on predictor type 									
plot_dt_primarycancer = plot_dt_select[type == "Primary cancer"]
plot_dt_normal_tissue = plot_dt_select[type == "Normal tissue"]									
plot_dt_cancer_cell_line = plot_dt_select[type == "Cancer cell line"]
									
plot_dt_RT = plot_dt_select[type == "RT"]
phases = unlist(lapply(as.character(plot_dt_RT$variable), 
					   function(x) unlist(strsplit(x, split = " "))[2]))
plot_dt_RT_early = plot_dt_RT[which(phases %in% c("G1", "S1", "S2"))]
plot_dt_RT_late = plot_dt_RT[which(phases %in% c("S3", "S4", "G2"))]

									
#Create faceted plot for different predictor types (data.table long-form)
plot_dt_facet = as.data.table(rbind.data.frame(plot_dt_primarycancer,
											   plot_dt_normal_tissue,
											   plot_dt_cancer_cell_line,
											   plot_dt_RT_early,
											   plot_dt_RT_late))
									
plot_dt_facet$facet = c(rep("Primary cancer", nrow(plot_dt_primarycancer)),
						rep("Normal tissue", nrow(plot_dt_normal_tissue)),
						rep("Cancer cell line", nrow(plot_dt_cancer_cell_line)),
						rep("RT Early", nrow(plot_dt_RT_early)), 
						rep("RT Late", nrow(plot_dt_RT_late)))

#Plot SHAP vs. CA density plot					   
plot_dt_facet$facet = factor(plot_dt_facet$facet, 
							 levels = c("Primary cancer", "Normal tissue", "RT Early", "RT Late"))									
pdf(pff("data/003H_SHAP_facet_densityplot.pdf"), height = 3, width = 8.5)                                      
ggplot(plot_dt_facet, aes(y = Input, x = SHAP))+
	facet_wrap(~facet, scales = "free", nrow = 1)+ #Facet based on predictor type
    geom_hex(bins = 50)+ #Add hexagonal desnity plot
    geom_smooth(method = "loess", span = 0.6, se = F)+ #Add loess smooth line
    theme_bw()+
    scale_fill_gsea()+ #Control colours
	labs(y = "Input(log1p, scaled)") #Add labels
dev.off()

#Correlations for the different plots
cor(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "spearman")
# [1] -0.7366388
cor.test(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "spearman")$p.val
# [1] 0
cor(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "spearman")
# [1] -0.7855442
cor.test(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "spearman")$p.val
# [1] 0	
cor(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "spearman")
# [1] -0.8375119
cor.test(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "spearman")$p.val
# [1] 0									
cor(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "spearman")
# [1] 0.7515782
cor.test(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "spearman")$p.val
# [1] 0													
									
									