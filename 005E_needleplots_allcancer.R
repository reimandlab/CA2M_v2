source("000_HEADER.R")

#Import liftover files for hg19 to hg38
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

#Get ENCODE ATAC-Seq bigwig paths
ENCODE_ATAC_track_names = colnames(fread(pff("data/001D_ENCODE_ATACSeq_100KBwindow_processed.csv"))[,-c(1,2)])
ENCODE_ATAC_track_paths = list.files(paste0(input_data_dir, "ENCODE_ATACSeq/"), 
						 full.names = T)[-c(1,2,208)]
#Get TCGA ATAC-Seq bigwig paths
TCGA_ATAC_track_paths = list.files("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/TCGA_ATACSEQ_Individual_Bigwigs/", 
								   full.names=T)
TCGA_ATAC_track_names =  gsub("_","-", unlist(lapply(list.files("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/TCGA_ATACSEQ_Individual_Bigwigs/"),
	function(x) unlist(strsplit(x,split = ".", fixed = T))[1])))
													 
#Get TCGA ATAC-Seq patients
matching_table = fread(paste0(input_data_dir, "TCGA_ATACSEQ_identifier_mapping.txt"))
ATAC_patients_all = unlist(lapply(matching_table$Case_ID, 
								 function(x) paste(unlist(strsplit(x, split = "-"))[1:3], collapse = "-")))
TCGA_ATAC_track_patients = ATAC_patients_all[match(TCGA_ATAC_track_names, matching_table$bam_prefix)]

#Get ENCODE Replication Timing paths							  
RT_paths = list.files("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2MutData/RepliSeq_Tracks/",full.names=T)[-c(97,98)]
RT_names = unlist(lapply(list.files("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2MutData/RepliSeq_Tracks/"),
	function(x) unlist(strsplit(x,split=".bigWig"))[1]))[-c(97,98)]
	
#Get GEO CLL ATAC-Seq bigwig paths
CLL_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/CLL_tracks/"), full.names = T)
CLL_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/CLL_tracks/"))

#Get GEO Prostate ATAc-Seq bigwig paths
Prostate_paths = list.files(paste0(input_data_dir, 
								   "GEO_ATACSeq_tracks/Prostate_tracks/"), 
							full.names = T)
Prostate_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Prostate_tracks/"))
						 
#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           
					 
#Load in PCAWG MAF file
PCAWG_mutations_dt_hg38 = fread(pff("data/001A_PCAWG_mutations_hg38.csv"))
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))[,-c(1, 2)][-c(2858:2866, 20254, 24212:24219)]
									  
#Load in window qvals from 100KB-scale analysis
q_val_dt = fread(pff("data/005C_qvaldt_100KBerrorwindows.csv"))

#Load in human genome hg38
genome = BSgenome.Hsapiens.UCSC.hg38
									  
#Load in predictor importances from 1MB-analysis for all cohort
importance = fread(pff("data/003A_feature_importances_allpreds.csv"))		

#Get paths for model errors from 100KB-scale analysis
error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
#Get bootsrap test results for 1MB-scale analysis
bootstrap_paths = list.files(pff("data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types = unlist(lapply(list.files(pff("data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
#Load in permutation test paths for 1MB-scale analysis
p_val_paths = list.files(pff("data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                           							  
#Function to get mutations and ATAC-Seq every 1kb for 100kb windows of interest
get_plot_data = function(cohort_name){
	
	
	
	print(cohort_name)
	
	#Get highest error (most significant) window for cohort 
	window_index = which.min(q_val_dt[[cohort_name]])
	
	#Get 100KB window of interest as GRanges
	gr.window = GRanges(q_val_dt$chr[window_index], 
						IRanges(q_val_dt$start[window_index], q_val_dt$start[window_index] + 99999))
	
	#Get 1kb windows for whole chromosome
    gr.windows_small = tileGenome(seqinfo(Hsapiens)[q_val_dt$chr[window_index]], 
							tilewidth = 1000, 
							cut.last.tile.in.chrom = TRUE)
	
	#Get 1kb windows overlapping 100KB window of interest
	gr.windows = gr.windows_small[subjectHits(findOverlaps(gr.window, gr.windows_small))]
	
	###Get mutations for each 100bp window
	
	#Load in mutation data
    mut = PCAWG_mutations_dt_hg38[Project_Code == cohort_name]
    
    #Convert mutations file into GRanges object
    mutation.gr = GRanges(mut$Chromosome, 
						  IRanges(mut$Start_position,mut$End_position))
    
    #Get number of counts of mutations overlapping each genomic 1kb window in 100KB window of interest
    mutation_counts = as.numeric(table(factor(queryHits(findOverlaps(gr.windows, mutation.gr)), 
											  levels = 1:length(gr.windows))))
	
	
	###Get ATAC-seq for each 100bp window for best predictors
	
	#Get predictors for cohort
	predictors = importance[[1]]

	#Load in predictor importances
    data = importance[[cohort_name]]
	
	#Get permutation test p-values for each predictor
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))
		
	#Annotate significant predictors								
    significance = ifelse(predictor_pvals == 0, "*", "")
    names(significance) = predictors
	
	#Get predictor importance means from bootstrap experiment for significant predictors
    Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
														   full.name = T), fread)))
    sig_predictor_means = unlist(lapply(which(significance == "*"), 
										function(x) mean(Bootstrap_dt[[x]])))
	
	#Get top 5 predictors based on importance mean from bootstrap test									
    top_preds = names(sig_predictor_means[order(sig_predictor_means, 
													 decreasing = T)][1:5])	
	#Get average CA/RT from predictor bigwig for each top predictor
	get_pred_data = function(top_pred){
		
		#If predictor is in TCGA get path and import bigwig
		if(top_pred %in% predictors[1:381]){
			
			#Get path
			path = TCGA_ATAC_track_paths[which(TCGA_ATAC_track_patients == unlist(strsplit(top_pred, split = " ", fixed = T))[2])][1]
			
			#Import bigwig
			BigWig = import(path)

		#If predictor is in ENCODE get path and import bigwig
		}else if(top_pred %in% predictors[382:584]){
			
			#Get path
			path = list.files(paste0(ENCODE_ATAC_track_paths[which(ENCODE_ATAC_track_names == top_pred)], 
									 "/released/GRCh38/fold_change_over_control/bigwig/rep1/"), full.names = T)[1]
			
			#Import bigwig
			BigWig = import(path)

		#If predictor is in RT get path, import bigwig, and lift over to hg38
		}else if(top_pred %in% predictors[772:867]){
			
			#Get path
			path = RT_paths[which(predictors[772:867] == top_pred)]
			
			#Import bigwig
			BigWig = import(path)

			#Liftover to hg38
			BigWig = unlist(liftOver(BigWig, ch))
			
			
		#If predictor is in GEO CLL predictors get path, import bigwig, and lift over to hg38	
		}else if(top_pred %in% predictors[700:754]){
			
			#Get path
			path = CLL_paths[which(CLL_names == top_pred)]
			
			#Import bigwig
			BigWig = import(path)

			#Liftover to hg38
			BigWig = unlist(liftOver(BigWig, ch))
			
		#If predictor is in GEO prostate predictors get path, import bigwig, and lift over to hg38	
		}else if(top_pred %in% predictors[762:771]){
			
			#Get path
			path = Prostate_paths[which(Prostate_names == top_pred)]
			
			#Import bigwig
			BigWig = import(path)

			#Liftover to hg38
			BigWig = unlist(liftOver(BigWig, ch))
		}

		#Get an RleList object with one coverage vector per seqlevel in the BigWig
		gr.data.new = GenomicRanges::coverage(BigWig, weight="score")

		#Match seqlevels between input and human genome
		seqlevels(gr.windows, pruning.mode="coarse") = names(gr.data.new)

		#Get binned average CA/RT for each window
		gr.data.new = binnedAverage(gr.windows, gr.data.new, "value")
		
		return(gr.data.new$value)}
	
	#Get average CA/RT for small windows for top predictors									
	pred_data = as.data.table(do.call("cbind.data.frame", lapply(top_preds, get_pred_data)))
	
	#Get top predictor descriptions
	top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_preds, predictor_descriptions$Predictor_names)]	
	colnames(pred_data) = top_predictor_descriptions
	
	#Get observed and expected mutations (and FDR) for 100KB window of interest
	error_dt = fread(error_dt_paths[which(error_dt_cancertypes == cohort_name)])
	obs = error_dt$observed[window_index]
	exp = error_dt$predicted[window_index]
	FDR = q_val_dt[[cohort_name]][window_index]
	

	#Combine results into one data frame
	data = cbind.data.frame(chr = as.character(seqnames(gr.windows)),
							start = start(gr.windows),
							mut = mutation_counts,
							pred_data,
							cohort_name = rep(cohort_name, length(mutation_counts)),
							obs = rep(obs, length(mutation_counts)),
							exp = rep(exp, length(mutation_counts)),
							FDR = rep(FDR, length(mutation_counts)))
	
	
	#Save results
	fwrite(data, pff(paste0("data/005H_needle_plot_data_ct/", cohort_name, ".csv")))
}

# cancer_types_to_check = colnames(q_val_dt)[-c(1,2)][-5][18:22]
										
#Get needleplot data for all cancer types of interest										
lapply(colnames(q_val_dt)[-c(1,2)], get_100bp_data)
