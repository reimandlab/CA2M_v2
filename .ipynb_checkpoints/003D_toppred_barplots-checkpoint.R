source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

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

plot_barplots = function(cohort_index, top_N_predictors){
	
	print(cohort_index)
	cohort_name = colnames(importance)[-1][cohort_index]
	predictors = importance[[1]]

	#Get permutation test p-values for each predictor
    data = importance[[cohort_index + 1]]
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))
									
    significance = ifelse(predictor_pvals == 0, "*", "")
    names(significance) = predictors
	top_N_predictors = min(top_N_predictors, sum(significance == "*"))
	
	#Get predictor importance means and sd from bootstrap experiment for significant predictor                     
    Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
														   full.name = T), fread)))
    sig_predictor_means = unlist(lapply(which(significance == "*"), 
										function(x) mean(Bootstrap_dt[[x]])))                                   
    top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
													 decreasing = T)][1:top_N_predictors])
    top_predictor_means = as.numeric(sig_predictor_means[order(sig_predictor_means, 
															   decreasing = T)][1:top_N_predictors])
    top_predictor_SDs = unlist(lapply(top_predictors, function(x) sd(Bootstrap_dt[[x]])))								

    #Get predictor categories
    top_predictor_categories = paste0(predictor_descriptions$Predictor_categories[match(top_predictors, 
																	   predictor_descriptions$Predictor_names)], " CA")
	top_predictor_categories = ifelse(top_predictor_categories=="RT CA", 
									  "Replication timing", top_predictor_categories)
									  
    top_predictor_matching = ifelse(top_predictors %in% matching_list[[cohort_name]], 
									"Matching", 
									"Non-Matching")
    top_predictor_fill = unlist(lapply(1:length(top_predictors), 
									   function(x) ifelse(top_predictor_matching[x]=="Matching", 
														  paste(top_predictor_categories[x], "from matching tissue", sep = " "), 
														  paste(top_predictor_categories[x], " from non-matching tissue", sep = ""))))	
									   
	top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_predictors, predictor_descriptions$Predictor_names)]	

	if(any(top_predictor_descriptions == "")){
		top_predictor_descriptions[which(top_predictor_descriptions == "")] = top_predictors[which(top_predictor_descriptions == "")]
	}								   
	#Create data table for plotting                         
	results.m = as.data.table(cbind.data.frame(value = top_predictor_means, 
										   variable = top_predictors, 
										   SD = top_predictor_SDs, 
										   significance = rep("*", length(top_predictors)), 
										   fill = top_predictor_fill))
    
    results.m$variable = factor(results.m$variable, levels = top_predictors)
    results.m$fill = factor(results.m$fill, levels = c("Primary cancer CA from matching tissue", 
													   "Normal tissue CA from matching tissue", 
													   "Cancer cell line CA from matching tissue",
													   "Replication timing from matching tissue",
													   "Primary cancer CA from non-matching tissue", 
													   "Normal tissue CA from non-matching tissue", 
													   "Cancer cell line CA from non-matching tissue",
													   "Replication timing from non-matching tissue"))    
    
    #Get max/min value and SD for plot dims
    max_bar = max(results.m$value)
    max_SD = results.m[value == max_bar]$SD
    
    min = min(results.m$value - results.m$SD)
    if(min >= 0){min = 0}								 
	
	#Set fill colors								 
	fill_colours = c("#FFC300", "#58D68D", "#2980B9", "#FF3352",
					 "#9A7d0A", "#186A3B", "#154360", "#641E16")
									   
					

	#Plot barplot without legend
	p = ggplot(results.m, aes(x = variable, 
							  y = value, 
							  label = significance, 
							  fill = fill))+
		geom_bar(stat = "identity", color = "black")+
		geom_text(aes(y=value+SD+SD/8), colour = "black", size = 6)+
		geom_errorbar(aes(ymin = value - SD, ymax = value + SD), width=.2,
				 position=position_dodge(.9))+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.5, "cm"))+
		labs(x = "Predictor", 
			 y = "Importance (Mean increase MSE)", 
			 fill = "Predictor category",
			 title = cohort_name)+
		scale_fill_manual(drop = F, values = fill_colours)+
		scale_x_discrete(labels = top_predictor_descriptions)
		
}

pdf(pff("/data/003D_top_pred_barplots.pdf"), width = 14)
lapply(c(1:26)[-8], plot_barplots, 5)
dev.off()