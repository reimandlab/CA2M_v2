source("000_HEADER.R")

#Load in sample size dataset
#Load in sample size dataset
PCAWG_mutations_dt_hg38_SNP = fread(paste0(data_dir, "data/001G_PCAWG_MAF_SNP_withSBSSIG.csv"))
PCAWG_sample_table = PCAWG_mutations_dt_hg38_SNP[, 
												 .(Cohort_size = uniqueN(Donor_ID)), 
												 by = Project_Code][order(Cohort_size, decreasing = T)]


#Load in mutation tracks 
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]
mut_burden = apply(mutation_counts_dt, 2, mean)
mut_burden = mut_burden / PCAWG_sample_table$Cohort_size[match(names(mut_burden), PCAWG_sample_table$Project_Code)]

#Load in model accuracies
model_results = fread(pff("data/003A_feature_importances_allpreds.csv"))
model_accuracies = model_results[1, -1]

cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")
									
#Plot cancer sample size vs model accuracy
plot_dt = as.data.table(cbind.data.frame(cohort_name = cancer_types_to_keep, 
										mut_burden = mut_burden[match(cancer_types_to_keep, names(mut_burden))], 
										model_accuracy = as.numeric(model_accuracies)[match(cancer_types_to_keep, colnames(model_accuracies))]))	
									
rho = cor(plot_dt[[2]], plot_dt[[3]], method = "spearman")
									
pdf(pff("data/SF5_mutburden_modelaccuracy_scatterplot.pdf"))
ggplot(plot_dt, aes(x = mut_burden, 
					y = model_accuracy, 
				    label = cohort_name))+
	geom_point()+
	theme_bw()+
	theme()+
	labs(x = "Cohort mutation burden/MB", y = "CA + RT model accuracy")+
	geom_text_repel()+
	annotate("text", label = paste0("r=", round(rho, 2)), color = "red", x = 1, y = 0.97, size = 5)+
	ylim(c(0,1))																	
dev.off()									