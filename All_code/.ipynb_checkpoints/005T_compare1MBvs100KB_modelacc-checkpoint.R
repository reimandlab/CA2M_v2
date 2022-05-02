source("000_HEADER.R")

#Load 1MB SBS binned mutation count dataset paths
SBS_1MB_cohorts = gsub(".csv", "", list.files("data/001I_PCAWG_sigs_new"))
SBS_1MB_paths = list.files("data/001I_PCAWG_sigs_new", full.names = T)

#Load 100KB binned mutation count dataset
mutation_counts_dt_100KB = fread(pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))[,-c(1, 2)] 

#Load 100KB SBS binned mutation count dataset paths
SBS_100KB_cohorts = gsub(".csv", "", list.files("data/005P_mutsigs_100KB"))
SBS_100KB_paths = list.files("data/005P_mutsigs_100KB", full.names = T)

#Core cancer types to keep
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", 
						 "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", 
						 "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	


#Load in SBS model accuracies at 1MB
SBS_1MB_accuracies_cohorts = gsub(".csv", "", list.files("data/004A_Sig_RF_Results"))
SBS_1MB_accuracies_paths = list.files("data/004A_Sig_RF_Results", full.names = T)

#Load in SBS model accuracies at 1MB
SBS_100KB_accuracies_cohorts = gsub(".csv", "", list.files("data/005S_Sig_RF_Results_100KB"))
SBS_100KB_accuracies_paths = list.files("data/005S_Sig_RF_Results_100KB", full.names = T)

#Load in all mutation model accuracy
model_accuracy_1MB = fread(pff("data/003A_feature_importances_allpreds.csv"))
model_accuracy_100KB_paths = list.files("data/005A_100KB_RF_errorDT_results", full.names = T)
model_accuracy_100KB_cohorts = gsub(".csv", "", list.files("data/005A_100KB_RF_errorDT_results"))

#Function to convert R2 to adjusted R2
adjust_R2 = function(R2, n, k){
	adjusted_R2 = 1 - (1 - R2)*(n - 1)/(n - k - 1)
	return(adjusted_R2)
}

#Function to get percent of windows with no mutations for 1MB dataset
get_model_accuracy = function(cohort, scale){
	
	if(scale == "1MB"){
		SBS_paths = SBS_1MB_paths
		SBS_cohorts = SBS_1MB_cohorts
		
		model_paths = SBS_1MB_accuracies_paths
		model_cohorts = SBS_1MB_accuracies_cohorts
		
		All_accuracy = model_accuracy_1MB[[cohort]][1]
	}else{
		SBS_paths = SBS_100KB_paths
		SBS_cohorts = SBS_100KB_cohorts
		
		model_paths = SBS_100KB_accuracies_paths
		model_cohorts = SBS_100KB_accuracies_cohorts
		
		All_accuracy = fread(model_accuracy_100KB_paths[which(model_accuracy_100KB_cohorts == cohort)])
		All_accuracy = cor(All_accuracy$observed, All_accuracy$predicted)**2
		All_accuracy = adjust_R2(All_accuracy, 24440, 869)
	}
	
	#Load in SBS mutation dataset
	SBS_muts = fread(SBS_paths[which(SBS_cohorts == cohort)])[,-c(1, 2)]
	
	#Keep only SBS with at least 20k mutations
	SBS_muts = SBS_muts[,.SD, .SDcols = which(colSums(SBS_muts)>20000)]
	
	#Load in model accuracies
	SBS_accuracies = fread(model_paths[which(model_cohorts == cohort)])[1, -1]
	
	#Keep only SBS with 20k mutations
	SBS_accuracies = SBS_accuracies[,.SD,.SDcols = colnames(SBS_muts)]

	#Create data frame of results
	results_df = cbind.data.frame(cohort = rep(cohort, nrow(SBS_accuracies)), 
								  SBS = colnames(SBS_accuracies), 
								  adj_R2 = as.numeric(SBS_accuracies))
	
	#Add ALL mutations from cohort to data frame
	results_df = rbind.data.frame(results_df, c(cohort, "ALL", All_accuracy))
							 
	return(results_df)}

Accuracies_1MB = as.data.table(do.call("rbind.data.frame", lapply(cancer_types_to_keep, get_model_accuracy, "1MB")))							 
Accuracies_100KB = as.data.table(do.call("rbind.data.frame", lapply(cancer_types_to_keep, get_model_accuracy, "100KB")))	


plot_dt = rbind.data.frame(Accuracies_1MB, Accuracies_100KB) 						 
plot_dt$scale = c(rep("1MB", nrow(Accuracies_1MB)), 
				  rep("100KB", nrow(Accuracies_100KB)))		

plot_dt$adj_R2 = as.numeric(plot_dt$adj_R2)
SBS_order = c("ALL", "SBS1", "SBS5", "SBS40", 
			  "SBS2", "SBS13", 
			  "SBS3", "SBS6","SBS26", "SBS36", 
			  "SBS4", "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS22","SBS29", "SBS35", 
			  "SBS9", "SBS12", "SBS16","SBS17a", "SBS17b", "SBS18","SBS41", "SBS44")			


plot_dt = plot_dt[SBS %in% SBS_order]						 
plot_dt$SBS = factor(plot_dt$SBS, levels = 	SBS_order)
plot_dt$scale = factor(plot_dt$scale, levels = 	c("1MB", "100KB"))
plot_dt$point = paste0(plot_dt$cohort, plot_dt$SBS)
plot_dt = plot_dt[point!= "Lung-AdenoCASBS3"] #Remove SBS in 100KB that's not in 1MB


#Create barplots
Barplot = ggplot(na.omit(plot_dt),
					 aes(x = SBS, 
						 y = adj_R2,
						 fill = scale))+
				facet_wrap(~cohort, ncol = 2)+
				geom_bar(stat='identity', position='dodge', color = "black")+
				theme_bw()+
				theme(axis.title = element_text(size=9), #Format text
					axis.text.x = element_text(size = 8, 
											   angle = 90, 
											   hjust = 1, 
											   vjust = 0.5,
											   colour = "black"),
					axis.text.y = element_text(size = 8, colour = "black"),
					 legend.text = element_text(size = 8),
					 legend.title =element_text(size = 10))+
				labs(x = "Cohort", fill = "Scale", y = "Model Accuracy (Adj. R2)")	 
							 
pdf(pff("data/005T_1MBvs100KB_accuracybarplot.pdf"), width = 12, height = 15)	
Barplot
dev.off()


#Paired boxplots
plot_dt_boxplot = plot_dt[SBS %in% c("ALL", "SBS1", "SBS5")]
Boxplot = ggplot(plot_dt_boxplot, aes(x = scale, y = adj_R2)) +
			facet_wrap(~SBS, nrow = 1)+
			geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70") +
			geom_point(colour="red", size=2, alpha=0.5) +
			geom_line(aes(group=cohort), colour="red", linetype="11") +
			theme_bw()+
			labs(x = "Scale" , y = "Adj. R2")
				 
pdf(pff("data/005T_1MBvs100KB_accuracyboxplot.pdf"))	
Boxplot
dev.off()	

#Run paired wilcoxon test
p_val_all = wilcox.test(adj_R2 ~ scale, data = plot_dt_boxplot[SBS == "ALL"], paired = TRUE, alternative = "g")$p.val
p_val_SBS1 = wilcox.test(adj_R2 ~ scale, data = plot_dt_boxplot[SBS == "SBS1"], paired = TRUE, alternative = "g")$p.val
p_val_SBS5 = wilcox.test(adj_R2 ~ scale, data = plot_dt_boxplot[SBS == "SBS5"], paired = TRUE, alternative = "g")$p.val


#Remove SBS found in 100KB and not in 1MB

Boxplot = ggplot(plot_dt, aes(x = scale, y = adj_R2)) +
			geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70") +
			geom_point(colour="red", size=2, alpha=0.5) +
			geom_line(aes(group=point), colour="red", linetype="11") +
			theme_bw()+
			labs(x = "Scale" , y = "Adj. R2")
				 
pdf(pff("data/005T_1MBvs100KB_accuracyboxplot_allSBS.pdf"))	
Boxplot
dev.off()	

p_val_all = wilcox.test(plot_dt[scale == "1MB"]$adj_R2, plot_dt[scale == "100KB"]$adj_R2, paired = TRUE, alternative = "g")$p.val


#Create final boxplots
plot_dt$facet = ifelse(plot_dt$SBS == "ALL", "ALL", "SBS")

#Get p-vals
p_val_all = wilcox.test(plot_dt[facet == "ALL" & scale == "1MB"]$adj_R2, plot_dt[facet == "ALL" & scale == "100KB"]$adj_R2, paired = TRUE, alternative = "g")$p.val

p_val_SBS = wilcox.test(plot_dt[facet == "SBS" & scale == "1MB"]$adj_R2, plot_dt[facet == "SBS" & scale == "100KB"]$adj_R2, paired = TRUE, alternative = "g")$p.val

p_val_dt = cbind.data.frame(facet = c("ALL", "SBS"),
							label = c(p_val_all, p_val_SBS))
							

Boxplot = ggplot(plot_dt, aes(x = scale, y = adj_R2)) +
			facet_wrap(~facet, nrow = 1)+
			geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70", notch = T) +
			geom_point(colour="red", size=2, alpha=0.5) +
			geom_line(aes(group=point), colour="red", linetype="11") +
			theme_bw()+
			labs(x = "Scale" , y = "Adj. R2")+
			geom_text(
			  data    = p_val_dt,
			  mapping = aes(x = -Inf, y = -Inf, label = paste0("p=",label)),
			  hjust   = -0.1,
			  vjust   = -1)

pdf(pff("data/005T_1MBvs100KB_accuracyboxplot_allvsSBS.pdf"))	
Boxplot
dev.off()	
