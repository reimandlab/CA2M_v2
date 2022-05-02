source("000_HEADER.R")

#Load 1 MB binned mutation count dataset
mutation_counts_dt_1MB = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

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


#Function to get percent of windows with no mutations for 1MB dataset
get_percent_no_mut = function(cohort, scale){
	
	if(scale == "1MB"){
		SBS_paths = SBS_1MB_paths
		SBS_cohort = SBS_1MB_cohorts
		Muts = mutation_counts_dt_1MB
	}else{
		SBS_paths = SBS_100KB_paths
		SBS_cohort = SBS_100KB_cohorts
		Muts = mutation_counts_dt_100KB}
	
	#Load in SBS mutation dataset
	SBS_muts = fread(SBS_paths[which(SBS_cohort == cohort)])[,-c(1, 2)]
	
	#Keep only SBS with at least 20k mutations
	SBS_muts = SBS_muts[,.SD, .SDcols = which(colSums(SBS_muts)>20000)]

	#Get percent of windows with no mutations
	Zero_mut_windows = apply(SBS_muts, 2, function(x) sum(x == 0)/length(x))
	
	#Create data frame of results
	results_df = cbind.data.frame(cohort = rep(cohort, length(Zero_mut_windows)), 
								  SBS = names(Zero_mut_windows), 
								  percent_nomut = as.numeric(Zero_mut_windows)*100)
	
	#Add ALL mutations from cohort to data frame
	results_df = rbind.data.frame(results_df, c(cohort, "ALL", sum(Muts[[cohort]] == 0)/length(Muts[[cohort]])))
							 
	return(results_df)}
							 
Percent_nomut_1MB = as.data.table(do.call("rbind.data.frame", lapply(cancer_types_to_keep, get_percent_no_mut, "1MB")))							 
Percent_nomut_100KB = as.data.table(do.call("rbind.data.frame", lapply(cancer_types_to_keep, get_percent_no_mut, "100KB")))	
							 
#Keep SBS of interest
SBS_order = c("ALL", "SBS1", "SBS5", "SBS40", 
			  "SBS2", "SBS13", 
			  "SBS3", "SBS6","SBS26", "SBS36", 
			  "SBS4", "SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS22","SBS29", "SBS35", 
			  "SBS9", "SBS12", "SBS16","SBS17a", "SBS17b", "SBS18","SBS41", "SBS44")			

Percent_nomut_1MB$SBS = factor(Percent_nomut_1MB$SBS, levels = rev(SBS_order))

Percent_nomut_100KB$SBS = factor(Percent_nomut_100KB$SBS, levels = rev(SBS_order))
							 							 						 
							 
#Create 1MB heatmap
Heatmap_1MB = ggplot(na.omit(Percent_nomut_1MB),
					 aes(x = cohort, 
						 y = SBS,
						 fill = as.numeric(percent_nomut)))+
				geom_tile()+
				scale_fill_gradient(low="yellow", high="red", limits = c(0, 50))+
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
				labs(x = "Cohort", y = "SBS", fill = "% of windows with no mutations",
					 title = "1MB-scale windows")
										 
#Create 100KB heatmap
Heatmap_100KB = ggplot(na.omit(Percent_nomut_100KB),
					 aes(x = cohort, 
						 y = SBS,
						 fill = as.numeric(percent_nomut)))+
				geom_tile()+
				scale_fill_gradient(low="yellow", high="red", limits = c(0, 50))+
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
				labs(x = "Cohort", y = "SBS", fill = "% of windows with no mutations",
					 title = "100KB-scale windows")
							 
							 
pdf(pff("data/005Q_1MBvs100KB_percentnomut_heatmap.pdf"))	
list(Heatmap_1MB, Heatmap_100KB)
dev.off()
							 
							 
#Create barplots
#Create 1MB barplots
							 
# Percent_nomut_1MB = Percent_nomut_1MB[SBS %in% c("ALL", "SBS1", "SBS5")]						 
# Percent_nomut_100KB = Percent_nomut_100KB[SBS %in% c("ALL", "SBS1", "SBS5")]						 
plot_dt = rbind.data.frame(Percent_nomut_1MB, Percent_nomut_100KB) 						 
plot_dt$scale = c(rep("1MB", nrow(Percent_nomut_1MB)), 
				  rep("100KB", nrow(Percent_nomut_100KB)))		
plot_dt = plot_dt[SBS %in% c("ALL", "SBS1", "SBS5")]						 
plot_dt$SBS = factor(plot_dt$SBS, levels = 	c("ALL", "SBS1", "SBS5"))
plot_dt$scale = factor(plot_dt$scale, levels = 	c("1MB", "100KB"))
							 
Barplot = ggplot(na.omit(plot_dt),
					 aes(x = cohort, 
						 y = as.numeric(percent_nomut),
						 fill = scale))+
				facet_wrap(~SBS, nrow = 3)+
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
				labs(x = "Cohort", fill = "Scale", y = "% of windows with no mutations")	 
							 
pdf(pff("data/005Q_1MBvs100KB_percentnomutbarplot.pdf"))	
Barplot
dev.off()