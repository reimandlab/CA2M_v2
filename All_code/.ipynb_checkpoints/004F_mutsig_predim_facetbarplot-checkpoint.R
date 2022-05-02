source("000_HEADER.R")

#Get predictor importance values for RF models predicting signature mutation counts
importance_paths = list.files(pff("data/004A_Sig_RF_Results"), 
					   full.names = T)
imp_cohort_names = unlist(lapply(list.files(pff("data/004A_Sig_RF_Results")),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))
						 
#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           

#Load predictor matching list to PCAWG cohorts
matching_list = readRDS(pff("data/001J_Matching_predictors_dictionary.RDS"))

#Load in predictor importance values for RF models predicting all mutations in cohort
importance_ALL = fread(pff("/data/003A_feature_importances_allpreds.csv"))
							 								 
#Load in permutation test paths for RF models predicting all mutations
p_val_paths_ALL = list.files(pff("/data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types_ALL = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))                     
#Load in bootstrap test paths for RF models predicting all mutations                                 
bootstrap_paths_ALL = list.files(pff("/data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types_ALL = unlist(lapply(list.files(pff("/data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
						 
#Get plot data for cohort						 
get_plot_data_ct = function(cohort_name, top_N_predictors){

	#Load in predictor importance values for RF models predicting signature-based mutation counts
	importances = fread(importance_paths[which(imp_cohort_names == cohort_name)])
    
	#Get predictor names
	predictors = importances[[1]][-1]
    
	#Get signatures
    Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), cohort_name, ".csv"))[,-c(1,2)]

	#Get signatures with > 10,000 mutations (and get data for models predicting all mutations as well)
	Sigs = c("ALL", colnames(Sigs[,.SD, .SDcols = which(colSums(Sigs)>10000)]))
	
	#Keep only signatures with > 20,000 mutations for plot (and get data for models predicting all mutations as well)
	sigs_to_run = c("ALL", colnames(Sigs[,.SD, .SDcols = which(colSums(Sigs)>20000)]))
	
	#Get plot data for signature
	get_plot_data_sig = function(signature){
		
		#Get plot data for ALL mutations
		if(signature == "ALL"){
			

			#Load in predictor importances for RF model predicting all mutations
			data = importance_ALL[[cohort_name]]
			
			#Load in permutation test predictor importance values for RF model predicting all mutations
			permutation_dt = as.data.table(do.call("rbind.data.frame", 
											 lapply(list.files(p_val_paths_ALL[which(p_val_cancer_types_ALL == cohort_name)], 
														   full.name = T), fread)))
			#Load in bootstrap test predictor importance values for RF model predicting all mutations
			bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
											 lapply(list.files(bootstrap_paths_ALL[which(bootstrap_cancer_types_ALL == cohort_name)], 
														   full.name = T), fread)))
		#Get plot data for specific signature
		}else{
			
			#Load in predictor importances for RF model predicting mutation counts from signature
			data = importances[[signature]][-1]

			#Get indices for bootstrap and permutation tests
			cohort_index = which(gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)] == cohort_name)
			sig_index = which(Sigs[-1] == signature)		
			
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
		}
				
		#Get p-values of predictors from permutation test for RF models predicting signature-based mutation counts
		predictor_pvals = unlist(lapply(1:length(predictors), 
										function(x)  sum(data[x] < permutation_dt[[x]])/1000))

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
		#Create data table for plotting                         
		results.m = as.data.table(cbind.data.frame(value = top_predictor_means, 
											   variable = top_predictors, 
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
							

		#Keep barplot facet order                                 
		results.m$cohort_name = rep(cohort_name, top_N_predictors)
		results.m$signature = rep(signature, top_N_predictors)								   
		results.m$order_infacet = match(results.m$variable, top_predictors)
		i <<- i + top_N_predictors
										   
		return(results.m)								   
		
	}

										   
	result_dt=as.data.table(do.call("rbind.data.frame", lapply(Sigs_to_keep, get_plot_data_sig)))
    return(result_dt)}

#Get data for core 17 cancer types
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")								   
i = 0
result_dt = as.data.table(do.call("rbind.data.frame", 
								lapply(cancer_types_to_keep, get_plot_data_ct, 5)))                                         
										   
#Order signatures							   
SBS_order = c("ALL", "SBS1", "SBS5", "SBS40", 
			  "SBS2", "SBS13", 
			  "SBS3", "SBS6","SBS26", "SBS36", 
			  "SBS4", "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS22","SBS29", "SBS35", 
			  "SBS9", "SBS12", "SBS16","SBS17a", "SBS17b", "SBS18","SBS41", "SBS44")			

#Keep relevant signatures for plot										   
plot_dt_select = result_dt[signature %in% SBS_order]	

#Order cohorts										   
plot_dt_select$cohort_name = factor(plot_dt_select$cohort_name, 
									levels = cancer_types_to_keep)
#Order signatures							   
plot_dt_select$signature = factor(plot_dt_select$signature,
								 levels = SBS_order)
										   
#Set fill colors								 
fill_colours = c("#FFC300", "#58D68D", "#2980B9", "#FF3352",
				 "#9A7d0A", "#186A3B", "#6C3483", "#922B21")
										   
#Plot faceted barplot
p = ggplot(plot_dt_select, aes(x = order_infacet, 
							   y = value, 
							   fill = fill))+
    geom_bar(stat = "identity", color="black")+
    facet_grid(vars(signature), vars(cohort_name), scales = "free")+                                     
    theme_bw()+
    theme(axis.text = element_blank(),
		  axis.ticks = element_blank(),
         strip.text.y.right = element_text(angle = 0, size = 6),
         strip.text.x.top = element_text(size = 5),
        legend.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 6, face = "bold"),
        legend.key.size = unit(0.4, "cm"))+
	scale_fill_manual(drop = F, values = fill_colours)+
    labs(x = "", y = "", fill = "Predictor category")                                     

#Save plot										   
pdf(pff("data/004F_mutsig_predimp_facetbarplot.pdf"), width = 15, height = 10)
p
dev.off()										   
										   
							 