source("000_HEADER.R")

#Load in PCAWG MAF file
PCAWG_mutations_dt_hg38_SNP = fread(pff("001A_PCAWG_mutations_hg38.csv"))

#Create lookup table for cohort sample sizes
PCAWG_sample_table = PCAWG_mutations_dt_hg38_SNP[, 
												 .(Cohort_size = uniqueN(Donor_ID)), 
												 by = Project_Code][order(Cohort_size, decreasing = T)]


#Load in CA + RT model accuracies
cohort_names = unlist(lapply(list.files(paste0(data_dir, "data/004A_Sig_RF_Results/")), 
										function(x) unlist(strsplit(x, 
																	split = ".csv", 
																	fixed = T))[1]))
data_paths_tumor_RT = list.files(paste0(data_dir, "data/004A_Sig_RF_Results"), full.names = T)
accuracies_tumor_RT = unlist(lapply(data_paths_tumor_RT, 
								 function(x) as.numeric(fread(x)[["ALL"]][1])))

									
#Plot cancer sample size vs model accuracy
plot_dt = as.data.table(cbind.data.frame(cohort_name = cohort_names, 
										sample_size = PCAWG_sample_table$Cohort_size[match(cohort_names, 
																						   PCAWG_sample_table$Project_Code)], 
										model_accuracy = accuracies_tumor_RT))	
plot_dt[cohort_name=="PANCAN"]$sample_size = nrow(PCAWG_mutations_dt_hg38_SNP)				
									
rho = cor(plot_dt[[2]][-21], plot_dt[[3]][-21], method = "spearman")
									
pdf(pff("data/002A_samplesize_modelaccuracy_scatterplot.pdf"))
ggplot(plot_dt[-21], aes(x = (sample_size), 
					y = model_accuracy, 
				    label = cohort_name))+
	geom_point()+
	theme_bw()+
	theme()+
	labs(x = "Cohort sample size", y = "CA + RT model accuracy")+
	geom_text_repel()+
	annotate("text", label = paste0("r=", round(rho, 2)), color = "red", x = 40, y = 0.97, size = 5)+
	ylim(c(0,1))																	
dev.off()									