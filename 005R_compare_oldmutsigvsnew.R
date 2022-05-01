source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load 1 MB binned mutation count dataset
mutation_counts_dt_1MB = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

#Load NEW SBS binned mutation count dataset paths
SBS_NEW_cohorts = gsub(".csv", "", list.files("data/001I_PCAWG_sigs_new"))
SBS_NEW_paths = list.files("data/001I_PCAWG_sigs_new", full.names = T)

#Load OLD SBS binned mutation count dataset paths
SBS_OLD_cohorts = list.files("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/210317/data/001I_Sig_RMVs")
SBS_OLD_paths = list.files("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/210317/data/001I_Sig_RMVs", full.names = T)

get_scatterplot_data = function(cohort){
	
	#Load in NEW SBS 
	SBS_NEW = fread(SBS_NEW_paths[which(SBS_NEW_cohorts == cohort)])[,-c(1,2)]
	
	#Remove undetected signatures (min 20k mutations in cohort)
	SBS_NEW = SBS_NEW[,.SD, .SDcols = which(colSums(SBS_NEW)>20000)]

	#Load in OLD SBS
	SBS_OLD = as.data.table(do.call("cbind.data.frame",
									lapply(list.files(SBS_OLD_paths[which(SBS_OLD_cohorts == cohort)], full.names = T), 
										   function(x) fread(x)[[3]])))
	colnames(SBS_OLD) = gsub(".csv", "", list.files(SBS_OLD_paths[which(SBS_OLD_cohorts == cohort)]))
	
	#Keep same SBS for both datasets
	SBS_overlapping = intersect(colnames(SBS_NEW), colnames(SBS_OLD))
	
	create_plot = function(SBS){
		
		plot_dt = cbind.data.frame(facet = rep(paste(cohort, SBS), nrow(SBS_OLD)), 
								   NEW = SBS_NEW[[SBS]], 
								   OLD = SBS_OLD[[SBS]])
		
		cor = round(cor(plot_dt$OLD, plot_dt$NEW, method = "spearman"), 2)
		cor_dt = cbind.data.frame(facet = paste(cohort, SBS), cor = cor)
		
		return(list(plot_dt, cor_dt))

	
	}
	
	plot_dt_cohort = as.data.table(do.call("rbind.data.frame", lapply(SBS_overlapping, 
																	  function(x) create_plot(x)[[1]])))
	cor_dt_cohort = as.data.table(do.call("rbind.data.frame", lapply(SBS_overlapping, 
														   function(x) create_plot(x)[[2]])))								   
	return(list(plot_dt_cohort, cor_dt_cohort))}																  
																	  
								
										   
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", 
						 "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", 
						 "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	

														   
														   
														   
plot_dt = as.data.table(do.call("rbind.data.frame", 
								lapply(cancer_types_to_keep, function(x) get_scatterplot_data(x)[[1]])))
										   
cor_dt = as.data.table(do.call("rbind.data.frame", 
								lapply(cancer_types_to_keep, function(x) get_scatterplot_data(x)[[2]])))										   
										   
p = ggplot(plot_dt, aes(x = OLD, 
						y = NEW))+
	rasterise(geom_point(fill  = "black"), dpi = 300)+
	facet_wrap(~facet, nrow = 10, scales = "free")+								   
									   
# 	annotate("text", x = 0, y = max(plot_dt$NEW), 
# 			 label = paste0("r=",as.character(cor)), color = "red")+
	theme_bw()+
	theme()+
	labs(x = "SBS with max probability", y = "Probablistic")+
	scale_x_continuous(limits = c(0, NA))+
	scale_y_continuous(limits = c(0, NA))+
	geom_text(
	  data    = cor_dt,
	  mapping = aes(x = -Inf, y = Inf, label = cor),
	  hjust   = -0.2,
	  vjust   = 1.2,
	  color = "red"
	)

									   
									   
pdf(pff("data/005R_SBS_oldvsnew_scatterplot.pdf"), width = 15, height = 15)										   
p									   									   
dev.off()										   
