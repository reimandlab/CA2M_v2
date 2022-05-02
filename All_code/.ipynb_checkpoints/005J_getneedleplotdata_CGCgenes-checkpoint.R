source("000_HEADER.R")
library(Gviz)

#Import liftover files
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

#Get ENCODE ATAC paths
ENCODE_ATAC_track_names = colnames(fread(pff("data/001D_ENCODE_ATACSeq_100KBwindow_processed.csv"))[,-c(1,2)])
ENCODE_ATAC_track_paths = list.files(paste0(input_data_dir, "ENCODE_ATACSeq/"), 
						 full.names = T)[-c(1,2,208)]
#Load in paths to ATAC_Seq bigwigs
TCGA_ATAC_track_paths = list.files("", 
								   full.names=T)
TCGA_ATAC_track_names =  gsub("_","-", unlist(lapply(list.files(""),
	function(x) unlist(strsplit(x,split = ".", fixed = T))[1])))
													 
#Get ATAC patients
matching_table = fread(paste0(input_data_dir, "TCGA_ATACSEQ_identifier_mapping.txt"))
ATAC_patients_all = unlist(lapply(matching_table$Case_ID, 
								 function(x) paste(unlist(strsplit(x, split = "-"))[1:3], collapse = "-")))
TCGA_ATAC_track_patients = ATAC_patients_all[match(TCGA_ATAC_track_names, matching_table$bam_prefix)]

#Get RT paths							  
RT_paths = list.files("",full.names=T)[-c(97,98)]
RT_names = unlist(lapply(list.files(""),
	function(x) unlist(strsplit(x,split=".bigWig"))[1]))[-c(97,98)]
	
#Get GEO CLL paths
CLL_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/CLL_tracks/"), full.names = T)
CLL_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/CLL_tracks/"))

#Get GEO CLL paths
Prostate_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Prostate_tracks/"), full.names = T)
Prostate_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Prostate_tracks/"))
						 
#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           
					 
#Load in MAF file
PCAWG_mutations_dt_hg38 = fread(pff("data/001A_PCAWG_mutations_hg38.csv"))
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))[,-c(1, 2)][-c(2858:2866, 20254, 24212:24219)]
									  
#Load in window qvals
q_val_dt = fread(pff("data/005C_qvaldt_100KBerrorwindows.csv"))

#Load in human genome hg38
genome = BSgenome.Hsapiens.UCSC.hg38
									  
#Cohort importances
importance = fread(pff("data/003A_feature_importances_allpreds.csv"))		
predictors = importance[[1]]

#Get error_dt paths
error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
#Get bootsrap results
bootstrap_paths = list.files(pff("data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types = unlist(lapply(list.files(pff("data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
#Load in permutation test paths
p_val_paths = list.files(pff("data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
								   
#Load in significant windows overlapping CGC genes
CGC_gene_windows.gr = readRDS(pff("data/005I_CGC_gene_window_gr.RDS"))
								   
#Get corresponding path of predictor
get_pred_data = function(top_pred, gr.windows){

	#If predictor is in TCGA
	if(top_pred %in% predictors[1:381]){
		path = TCGA_ATAC_track_paths[which(TCGA_ATAC_track_patients == unlist(strsplit(top_pred, 
																					   split = " ", 
																					   fixed = T))[2])][1]
		#Import bigwig
		BigWig = import(path)

	#If predictor is in ENCODE
	}else if(top_pred %in% predictors[382:584]){
		path = list.files(paste0(ENCODE_ATAC_track_paths[which(ENCODE_ATAC_track_names == top_pred)], 
								 "/released/GRCh38/fold_change_over_control/bigwig/rep1/"), 
						  full.names = T)[1]
		#Import bigwig
		BigWig = import(path)

	#If predictor is in RT
	}else if(top_pred %in% predictors[772:867]){

		path = RT_paths[which(predictors[772:867] == top_pred)]

		BigWig = import(path)

		#Liftover to hg38
		BigWig = unlist(liftOver(BigWig, ch))
	}else if(top_pred %in% predictors[700:754]){

		path = CLL_paths[which(CLL_names == top_pred)]

		BigWig = import(path)

		#Liftover to hg38
		BigWig = unlist(liftOver(BigWig, ch))

	}else if(top_pred %in% predictors[762:771]){

		path = Prostate_paths[which(Prostate_names == top_pred)]

		BigWig = import(path)

		#Liftover to hg38
		BigWig = unlist(liftOver(BigWig, ch))
	}
	#Get an RleList object with one coverage vector per seqlevel in the BigWig
	gr.data.new = GenomicRanges::coverage(BigWig, weight="score")

	#Match seqlevels between input and human genome
	seqlevels(gr.windows, pruning.mode="coarse") = names(gr.data.new)

	#Get binned average for each window
	gr.data.new = binnedAverage(gr.windows, gr.data.new, "value")

	return(gr.data.new$value)}	
								   
								   
save_data = function(CGC_gene_window_num){
	print(CGC_gene_window_num)
	window = CGC_gene_windows.gr[CGC_gene_window_num]
	
	window_index = which(q_val_dt$chr == seqnames(window) & q_val_dt$start == start(window))
	
	#Get 100bp windows inside 100kb window 
	gr.window = GRanges(q_val_dt$chr[window_index], 
						IRanges(q_val_dt$start[window_index], q_val_dt$start[window_index] + 99999))

    gr.windows_small = tileGenome(seqinfo(Hsapiens)[q_val_dt$chr[window_index]], 
							tilewidth = 1000, 
							cut.last.tile.in.chrom=TRUE)
	
	gr.windows = gr.windows_small[subjectHits(findOverlaps(gr.window,gr.windows_small))]
	
	#Get cancer types for which this window is significant
	cancer_types_sig = colnames(q_val_dt)[-seq(2)][which(q_val_dt[,-seq(2)][window_index] < 0.05)]
	
	#Remove piloastrocytome if in cancer types
	if("CNS-PiloAstro" %in% cancer_types_sig){
		cancer_types_sig = cancer_types_sig[-which(cancer_types_sig == "CNS-PiloAstro")]
	}
	save_data_cohort = function(cohort_name){
	print(cohort_name)
		#Load in mutation data
		mut = PCAWG_mutations_dt_hg38[Project_Code == cohort_name]

		#Convert mutations file into GRanges object
		mutation.gr = GRanges(mut$Chromosome, 
							  IRanges(mut$Start_position,mut$End_position))

		#Get number of counts of mutations overlapping each genomic window
		mutation_counts = as.numeric(table(factor(queryHits(findOverlaps(gr.windows,mutation.gr)), 
												  levels = 1:length(gr.windows))))


		#Get permutation test p-values for each predictor
		data = importance[[cohort_name]]
		permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
														   full.name = T), fread)))
		predictor_pvals = unlist(lapply(1:length(predictors), 
										function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))

		significance = ifelse(predictor_pvals == 0, "*", "")
		names(significance) = predictors

		#Get predictor importance means and sd from bootstrap experiment for significant predictor                     
		Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
											 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
															   full.name = T), fread)))
		sig_predictor_means = unlist(lapply(which(significance == "*"), 
											function(x) mean(Bootstrap_dt[[x]])))                                   
		top_preds = names(sig_predictor_means[order(sig_predictor_means, 
														 decreasing = T)][1:5])	
		
		
		pred_data = as.data.table(do.call("cbind.data.frame", lapply(top_preds, get_pred_data, gr.windows)))

		#Get top pred descriptions
		top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_preds, predictor_descriptions$Predictor_names)]	
		colnames(pred_data) = top_predictor_descriptions

		#Get obs/expected
		error_dt = fread(error_dt_paths[which(error_dt_cancertypes == cohort_name)])
		obs = error_dt$observed[window_index]
		exp = error_dt$predicted[window_index]
		FDR = q_val_dt[[cohort_name]][window_index]


		#Get data frame of results
		data = cbind.data.frame(chr = as.character(seqnames(gr.windows)),
									start = start(gr.windows),
									mut = mutation_counts,
									pred_data,
									cohort_name = rep(cohort_name, length(mutation_counts)),
									gene = window$gene,
									obs = rep(obs, length(mutation_counts)),
									exp = rep(exp, length(mutation_counts)),
									FDR = rep(FDR, length(mutation_counts)))


		fwrite(data, pff(paste0("data/005J_CGC_needle_plot_data/", 
								cohort_name, CGC_gene_window_num, 
								".csv")))
		}

	lapply(cancer_types_sig, save_data_cohort)
	
}								   
lapply(48:length(CGC_gene_windows.gr), save_data)                           		