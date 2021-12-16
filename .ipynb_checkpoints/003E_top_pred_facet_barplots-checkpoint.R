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

get_plot_data = function(cohort_index, top_N_predictors){
	
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
			
	if(length(unique(top_predictors)!=top_N_predictors)){
		top_predictors[duplicated(top_predictors)] = paste0(top_predictors[duplicated(top_predictors)],"2")
	}

	if(any(top_predictor_descriptions == "")){
		top_predictor_descriptions[which(top_predictor_descriptions == "")] = top_predictors[which(top_predictor_descriptions == "")]
	}	

	#Create data table for plotting                         
	results.m = as.data.table(cbind.data.frame(value = top_predictor_means, 
										   variable = top_predictors, 
										   SD = top_predictor_SDs, 
										   significance = rep("*", length(top_predictors)),
										   description = top_predictor_descriptions,	   
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
	
	#Keep barplot facet order                                 
    results.m$facet=rep(cohort_name, top_N_predictors)
    results.m$order=match(results.m$variable, top_predictors)+i
    results.m$order_infacet=match(results.m$variable, top_predictors)
    i <<- i+5
									   
	return(results.m)
}

i = 0
plot_dt = as.data.table(do.call("rbind.data.frame",
	lapply(seq(26)[-8],  get_plot_data, 5)))         									   

cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")									   
plot_dt_select = plot_dt[facet %in% cancer_types_to_keep]
plot_dt_select$facet = factor(plot_dt_select$facet, levels = cancer_types_to_keep)									   
									   
#Set fill colors								 
fill_colours = c("#FFC300", "#58D68D", "#2980B9", "#FF3352",
				 "#9A7d0A", "#186A3B", "#6C3483", "#922B21")

#Get x axis labels
labels = as.character(plot_dt_select$description[match(unique(plot_dt_select$order), 
													   plot_dt_select$order)])									   
names(labels) = unique(plot_dt_select$order)  									   
									   
#Plot barplot without legend
pdf(pff("/data/003E_top_pred_facet_barplots.pdf"), width = 12, height = 8)
ggplot(plot_dt_select, aes(x = factor(order), 
						  y = value, 
						  label = significance, 
						  fill = fill,
						  group = factor(order)))+
	facet_wrap(~facet, nrow = 3, scales = "free")+
	geom_bar(stat = "identity", color = "black")+
	geom_text(aes(y=value+SD+SD/8), colour = "black", size = 6)+
	geom_errorbar(aes(ymin = value - SD, ymax = value + SD), width=.2,
			 position=position_dodge(.9))+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=  8, colour = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 6, colour = "black"),
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        legend.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.8,0.1))+	
	labs(x = "Predictor", 
		 y = "Importance (Mean increase MSE)", 
		 fill = "Predictor category")+
	scale_fill_manual(drop = F, values = fill_colours, guide=FALSE)+
    scale_x_discrete(breaks = names(labels), labels = labels)+
	geom_blank(aes(y = value+SD+SD/2))
dev.off()
									   
#Plot barplot without legend									   
p1_legend = ggplot(plot_dt_select, aes(x = factor(order), 
						  y = value, 
						  label = significance, 
						  fill = fill,
						  group = factor(order)))+
	facet_wrap(~facet, nrow = 3, scales = "free")+
	geom_bar(stat = "identity", color = "black")+
	geom_text(aes(y=value+SD+SD/8), colour = "black", size = 6)+
	geom_errorbar(aes(ymin = value - SD, ymax = value + SD), width=.2,
			 position=position_dodge(.9))+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=  8, colour = "black"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 6, colour = "black"),
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        legend.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.5, "cm"))+	
		  labs(x = "Predictor", 
		 y = "Importance (Mean increase MSE)", 
		 fill = "Predictor category")+
	scale_fill_manual(drop = F, values = fill_colours)+
    scale_x_discrete(breaks = names(labels), labels = labels)+
	geom_blank(aes(y = value+SD+SD/2))
												   

legend_full <- cowplot::get_legend(p1_legend)

pdf(pff("/data/003E_top_pred_facet_barplots_legend.pdf"), width = 4, height = 4)
grid.newpage()             
grid.draw(legend_full)                                      
dev.off()              											   