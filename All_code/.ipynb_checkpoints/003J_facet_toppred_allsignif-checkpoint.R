source("000_HEADER.R")

#Load in predictor importances for normal RF models
importance = fread(pff("data/003A_feature_importances_allpreds.csv"))

#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           

#Load predictor matching list to PCAWG cohorts
matching_list = readRDS(pff("data/001J_Matching_predictors_dictionary.RDS"))

#Load in permutation test result paths and names
p_val_paths = list.files(pff("data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 
#Load in bootstrap test result paths and names
bootstrap_paths = list.files(pff("data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types = unlist(lapply(list.files(pff("data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
									   
#Function to plot barplots of importances for top N predictors for specific cohort
get_plot_data = function(cohort_name){

	#Get predictor importances from permutation test
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
	#Get predictor names
	predictors = colnames(permutation_importances_dt)
	
	#Get predictor importances
    data = importance[[cohort_name]][match(predictors, importance[[1]])]
	
	#Get permutation test p-values for each predictor (fraction of times permuted importance is more significant than actual importance)
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))							
									
	#Predictor are significant if they are more important all 1000 permuted values
    significance = ifelse(predictor_pvals == 0, "*", "")
    names(significance) = predictors
	
	#Get predictor importances from bootstrap test
    Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
														   full.name = T), fread)))
	#Get the importance mean from the boostrap test for all significant predictors								
    sig_predictor_means = unlist(lapply(which(significance == "*"), 
										function(x) mean(Bootstrap_dt[[x]])))                                   
    
	#Rank significant predictors by bootstrap mean
	top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
													 decreasing = T)])
										
	#Get predictor importance means and sd from bootstrap experiment for each significant predictor                     
    top_predictor_means = as.numeric(sig_predictor_means[order(sig_predictor_means, 
															   decreasing = T)])
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
	
	#Some predictor names are duplicated so just make these unique								   
	if(length(unique(top_predictors))!=length(top_predictors)){
		top_predictors[duplicated(top_predictors)] = paste0(top_predictors[duplicated(top_predictors)],"2")
	}
	
	#Some predictors have no descriptions so just keep original predictor name for these								   
	if(any(top_predictor_descriptions == "")){
		top_predictor_descriptions[which(top_predictor_descriptions == "")] = top_predictors[which(top_predictor_descriptions == "")]
	}	

	#Create data table for plotting                         
	results.m = as.data.table(cbind.data.frame(value = top_predictor_means, 
										   variable = top_predictors,
											description = top_predictor_descriptions,
										   SD = top_predictor_SDs, 
										   significance = rep("*", length(top_predictors)), 
										   fill = top_predictor_fill))
    
	#Order predictors
    results.m$variable = factor(results.m$variable, levels = top_predictors)
									   
	#Order fill levels								   
    results.m$fill = factor(results.m$fill, levels = c("Primary cancer CA from matching tissue", 
													   "Normal tissue CA from matching tissue", 
													   "Cancer cell line CA from matching tissue",
													   "Replication timing from matching tissue",
													   "Primary cancer CA from non-matching tissue", 
													   "Normal tissue CA from non-matching tissue", 
													   "Cancer cell line CA from non-matching tissue",
													   "Replication timing from non-matching tissue"))    
    
    
    #Get max/min value and SD for plot dimensions
    max_bar = max(results.m$value)
    max_SD = results.m[value == max_bar]$SD
    
    min = min(results.m$value - results.m$SD)
    if(min >= 0){min = 0}								 
	
	#Keep barplot facet order                                 
    results.m$facet=rep(cohort_name, length(top_predictors))
    results.m$order=match(results.m$variable, top_predictors)+i
    results.m$order_infacet=match(results.m$variable, top_predictors)
    i <<- i+length(top_predictors)
									   
	return(results.m)
}

i = 0 #Index for facet order
				
#Keep only core 17 cancer types									   
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
#Get plot data for all cancer types									   
plot_dt = as.data.table(do.call("rbind.data.frame",
	lapply(cancer_types_to_keep,  get_plot_data)))         									   

								   
plot_dt$facet = factor(plot_dt$facet, levels = cancer_types_to_keep)									   
									   
#Set fill colors								 
fill_colours = c("#FFC300", "#58D68D", "#2980B9", "#FF3352",
				 "#9A7d0A", "#186A3B", "#6C3483", "#922B21")

#Get x axis labels
labels = as.character(plot_dt$description[match(unique(plot_dt$order), 
													   plot_dt$order)])									   
names(labels) = unique(plot_dt$order)  									   
									   
#Plot barplot without legend
pdf(pff("/data/003J_top_pred_facet_barplots_allsignificant.pdf"), width = 12, height = 10)
ggplot(plot_dt, aes(x = factor(order), 
						  y = value, 
						  label = significance, 
						  fill = fill,
						  group = factor(order)))+
	facet_wrap(~facet, nrow = 4, scales = "free")+ #Facet wrap cancer types
	geom_bar(stat = "identity", color = "black")+ #Add barplots
	geom_text(aes(y=value+SD+SD/8), colour = "black", size = 6)+ #Add significance
	geom_errorbar(aes(ymin = value - SD, ymax = value + SD), width=.2, #Add errorbars
			 position=position_dodge(.9))+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=  8, colour = "black"), #Control text
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 6, colour = "black"),
        strip.text = element_text(size = 8, colour = "black", face = "bold"),
        legend.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.8,0.1))+	
	labs(x = "Predictor",  #Add labels
		 y = "Importance (Mean increase MSE)", 
		 fill = "Predictor category")+
	scale_fill_manual(drop = F, values = fill_colours, guide=FALSE)+ #Add fill colors and legend
    scale_x_discrete(breaks = names(labels), labels = labels)+ #Add x axis labels
	geom_blank(aes(y = value+SD+SD/2)) #Control y-axis scale
dev.off()
									   
#Calculate enrichment of primary cancers among the significant predictors

#Get names of all predictors									   
all_predictors = predictor_descriptions$Predictor_names
									   
#Get significant predictors in core 17 cancer types (some repeat)									   
significant_predictors = as.character(plot_dt$variable)
									   
#Get which predictors are primary cancer									   
cancer_CA_predictors = predictor_descriptions[Predictor_categories == "Primary cancer"]$Predictor_names

#Get which of significant predictors are from primary cancer samples									   
top_cancer = sum(significant_predictors %in% cancer_CA_predictors)
									   
#Get which of significant predictors are not from primary cancer samples									   									   
top_noncancer = length(significant_predictors) - top_cancer

#Get which of non-significant predictors are from primary cancer samples
nontop_cancer = length(cancer_CA_predictors) - top_cancer
									   
#Get which of non-significant predictors are not from primary cancer samples									   
nontop_noncancer = length(all_predictors) - top_cancer - top_noncancer - nontop_cancer									   									   

#Create contingency table
cont_table = matrix(c(top_cancer, top_noncancer, 
					  nontop_cancer, nontop_noncancer),
					nrow = 2)
rownames(cont_table) = c("Cancer", "Non-cancer")
colnames(cont_table) = c("Significant", "Non-significant")

#Get Fisher's exact test p-value and expected value for enrichment of primary cancer predictors 
# among top predictors in core 17 cancer types									   
p_val = fisher.test(cont_table,
					alternative = "g")$p.val
p_val												 
# [1] 8.681534e-08		
													 
Expected_value = (length(significant_predictors) * length(cancer_CA_predictors))/ length(all_predictors)
Expected_value	#Expected number of primary cancer predictors in significant

# 80.42117
									   
Observed_value = top_cancer
Observed_value #Observed number of primary cancer predictors in top 85

#111									   