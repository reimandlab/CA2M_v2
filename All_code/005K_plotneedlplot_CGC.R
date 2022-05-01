source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
library(Gviz)
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Load in CGC genes for each CGC gene window
CGC_gene_windows.gr = readRDS(pff("data/005I_CGC_gene_window_gr.RDS"))
CGC_genes = CGC_gene_windows.gr$gene

plot_data = function(data_path_index){
	
	#Load in data
	data = fread(paths[data_path_index])
	cohort_name = data$cohort_name[1]
	
	#Get CGC gene
	path_name = path_names[data_path_index]
	CGC_gene = CGC_genes[as.numeric(gsub(cohort_name, "", path_name))]	
	
	top_pred = data$pred[1]
	exp = round(data$exp[1], 2)
	obs = data$obs[1]
	FDR = signif(data$FDR[1], 2)
	
	gen = "hg38"
	chr = data$chr[1]
	itrack = IdeogramTrack(genome = gen, chromosome = chr)
	gtrack = GenomeAxisTrack()

	grtrack = GeneRegionTrack(paste0(input_data_dir, "gencode.v38.annotation.gtf"),
							  stacking = 'squish', stackHeight = 0.3, name = 'Genes')
	
	muts = t(data.matrix(cbind(data$mut, rep(exp/length(data$mut), length(data$mut)))))
	groups_mut = factor(c("Observed", "Expected"))
	
	dtrack1 = DataTrack(data = muts, groups = groups_mut, start = data$start,
						end = data$start+999, chromosome = chr, genome = gen, 
						name = paste0("Mutations/KB\n", cohort_name), type = c("histogram"), 
						cex.axis = 0.7)
	
	ATAC = data.matrix(data[,c(4:8)])
	#Normalize ATAC
	ATAC = t(apply(ATAC, 2, function(x) (x-min(x))/(max(x)-min(x))))+0.01
	groups_ATAC = factor(make.unique(colnames(data)[c(4:8)]), levels = make.unique(colnames(data)[c(4:8)]))
	
	dtrack2 = DataTrack(data = ATAC, groups = groups_ATAC, start = data$start,
						end = data$start+999, chromosome = chr, genome = gen, 
						name = "CA or RT/KB", type = c("smooth"), 
						span = 0.2)
	
	plotTracks(list(itrack, gtrack, grtrack, dtrack1, dtrack2), 
			   from =  data$start[1], to = data$start[nrow(data)],  
			   transcriptAnnotation = 'symbol', collapseTranscripts = 'meta', 
			   main = paste0(cohort_name, "\nObs=", obs, " Exp=", exp, " FDR=", FDR, "\n", CGC_gene), 
			   background.title = "darkblue", cex.main = 0.9, extend.right = 20000, 
			   extend.left = 20000)
	
}

paths = list.files(pff(paste0("data/005J_CGC_needle_plot_data/")), full.names = T)
path_names = gsub(".csv", "", list.files(pff(paste0("data/005J_CGC_needle_plot_data/"))))
				   
				   
pdf(pff("data/005K_CGCneedleplots.pdf"), width  = 5, height = 6)
lapply(seq_along(paths), plot_data)
dev.off()