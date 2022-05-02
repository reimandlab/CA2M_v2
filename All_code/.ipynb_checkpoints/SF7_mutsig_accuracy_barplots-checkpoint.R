source("000_HEADER.R")

importance_paths = list.files(pff("data/004A_Sig_RF_Results"), 
					   full.names = T)
imp_cohort_names = unlist(lapply(list.files(pff("data/004A_Sig_RF_Results")),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))
model_results = fread(pff("data/003A_feature_importances_allpreds.csv"))
model_accuracies = model_results[1, -1]
								 
get_plot_data = function(cohort_name){
    
	Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), cohort_name, ".csv"))[,-c(1,2)]
	Sigs_to_keep = colnames(Sigs[,.SD, .SDcols = which(colSums(Sigs)>20000)])

	importances = fread(importance_paths[which(imp_cohort_names == cohort_name)]) 
    signatures = colnames(importances)[-1]
    
    Adj_R2 = as.numeric(importances[1, -1])
	All_adj_R2 = as.numeric(model_accuracies)[which(names(model_accuracies) == cohort_name)]	
	plot_dt = as.data.table(cbind.data.frame(cohort_name = rep(cohort_name, length(signatures)), 
							   signatures, 
							   Adj_R2, 
							   ALL = rep(All_adj_R2, length(signatures))))
	
	plot_dt = plot_dt[signatures %in% Sigs_to_keep]
	plot_dt$order = match(plot_dt$Adj_R2, sort(plot_dt$Adj_R2, decreasing = T))+i
	i<<- i+length(signatures)
	return(plot_dt)}
							 
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")								   
i = 0
result_dt = as.data.table(do.call("rbind.data.frame", 
								lapply(cancer_types_to_keep, get_plot_data)))    
							 
etiology_dt = as.data.table(cbind.data.frame(
    etiology = c(rep("APOBEC enzyme activity", 2), 
				 rep("Defective DNA repair", 4), 
				 rep("SBS1", 1), 
				 rep("Exogenous/Carcinogen", 8),
				 rep("SBS5, SBS40", 2),
				 rep("Unknown/Other", 8)), 
    signature = c(c("SBS13", "SBS2"), 
				  c("SBS3", "SBS6", "SBS26", "SBS44"), 
				  c("SBS1"), 
				  c("SBS4", "SBS7a", "SBS7b", "SBS7c","SBS7d","SBS29", "SBS22", "SBS35"),
				  c("SBS5", "SBS40"),
				  c( "SBS9", "SBS12","SBS16", "SBS17a", "SBS17b", "SBS18", "SBS37", "SBS41"))))							 
                            
result_dt$etiologies = etiology_dt$etiology[match(result_dt$signature, 
											  etiology_dt$signature)]
                            
#Create boxplots
result_dt$etiologies = factor(result_dt$etiologies, 
								 levels = c("SBS1", 
											"APOBEC enzyme activity", 
											"Defective DNA repair", 
											"Exogenous/Carcinogen",
											"SBS5, SBS40",
											"Unknown/Other"))  
                            
labels = as.character(result_dt$signatures[match(sort(result_dt$order), 
														 result_dt$order)])
names(labels) = sort(result_dt$order)    
							 
p = ggplot(result_dt, aes(x = factor(order), 
					  y = Adj_R2, 
					  fill = etiologies))+
	facet_wrap(~factor(cohort_name, levels = cancer_types_to_keep), 
			   scales = "free_x", ncol = 3)+
	geom_bar(stat = "identity")+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
		  legend.position = "bottom")+					 
	scale_fill_d3(drop = F)+
	scale_x_discrete(breaks = names(labels), labels = labels)+
	geom_hline(aes(yintercept = ALL), colour = "red")+
	labs(x = "SBS Signature", y = "Model accuracy (adj. R2)", fill = "Signature etiology")+
	ylim(c(0,1))						 



							 
pdf(pff("data/SF7_mutsig_accuracy_barplots.pdf"), width = 8, height = 10)
p
dev.off()							 
							 