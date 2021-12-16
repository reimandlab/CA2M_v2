date_tag = "210426"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/210317/"

#Load in sample size dataset
PCAWG_mutations_dt_hg38_SNP = fread(paste0(data_dir, "data/001G_PCAWG_MAF_SNP_withSBSSIG.csv"))
PCAWG_sample_table = PCAWG_mutations_dt_hg38_SNP[, 
												 .(Cohort_size = uniqueN(Donor_ID)), 
												 by = Project_Code][order(Cohort_size, decreasing = T)]


#Load in mutation tracks 
mut_tracks = fread(paste0(data_dir, "data/001D_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]
mut_burden = apply(mut_tracks, 2, mean)
mut_burden = mut_burden / PCAWG_sample_table$Cohort_size[match(names(mut_burden), PCAWG_sample_table$Project_Code)]
mut_burden[which(names(mut_burden) == "PANCAN")] = mean(mut_tracks[["PANCAN"]]/sum(PCAWG_sample_table$Cohort_size))

#Load in CA + RT model accuracies
cohort_names = unlist(lapply(list.files(paste0(data_dir, "data/004A_Sig_RF_Results/")), 
										function(x) unlist(strsplit(x, 
																	split = ".csv", 
																	fixed = T))[1]))
data_paths = list.files(paste0(data_dir, "data/004A_Sig_RF_Results"), full.names = T)
accuracies = unlist(lapply(data_paths, 
								 function(x) as.numeric(fread(x)[["ALL"]][1])))

									
#Plot cancer sample size vs model accuracy
plot_dt = as.data.table(cbind.data.frame(cohort_name = cohort_names, 
										mutation_burden = as.numeric(mut_burden)[match(cohort_names, 
																						   names(mut_burden))], 
										model_accuracy = accuracies))	
									
rho = cor(plot_dt[[2]], plot_dt[[3]], method = "spearman")
									
pdf(pff("data/002B_mutationburden_modelaccuracy_scatterplot.pdf"))
ggplot(plot_dt, aes(x = (mutation_burden), 
					y = model_accuracy, 
				    label = cohort_name))+
	geom_point()+
	theme_bw()+
	theme()+
	labs(x = "Average mutations per sample in 1 Mbp window", y = "CA + RT model accuracy")+
	geom_text_repel()+
	annotate("text", label = paste0("r=", round(rho, 2)), color = "red", x = 1, y = 0.97, size = 5)+
	ylim(c(0,1))								
dev.off()									