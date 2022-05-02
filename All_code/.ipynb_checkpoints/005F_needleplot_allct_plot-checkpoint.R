source("000_HEADER.R")

#Create needleplot for windows of interest
plot_data = function(data_path){
	
	#Load in plot data
	data = fread(data_path)
	
	#Get cohort name
	cohort_name = data$cohort_name[1]
	
	#Get name of top predictor
	top_pred = data$pred[1]
	
	#Get and round expected mutations for window of interest
	exp = round(data$exp[1], 2)
	
	#Get observed mutations for window of interest
	obs = data$obs[1]
	
	#Get FDR at N significant digits
	FDR = signif(data$FDR[1], 2)
	
	#Define genome as hg38
	gen = "hg38"
	
	#Get chromosome of window
	chr = data$chr[1]
	
	#Create ideogram track for needle plot
	itrack = IdeogramTrack(genome = gen, chromosome = chr)
	
	#Create genome axis track for needle plot
	gtrack = GenomeAxisTrack()
	
	#Create gene region track for needle plot
	grtrack = GeneRegionTrack(paste0(input_data_dir, "gencode.v38.annotation.gtf"),
							  stacking = 'squish', stackHeight = 0.3, name = 'Genes')
	
	#Get mutation data table with observed and average expected mutations for window
	muts = t(data.matrix(cbind(data$mut, rep(exp/length(data$mut), length(data$mut)))))
	groups_mut = factor(c("Observed", "Expected"))
	
	#Create data track with observed and expected mutations
	dtrack1 = DataTrack(data = muts, groups = groups_mut, start = data$start,
						end = data$start+999, chromosome = chr, genome = gen, 
						name = paste0("Mutations/KB\n", cohort_name), type = c("histogram"), 
						cex.axis = 0.7)
	#Get CA/RT values
	CA_RT = data.matrix(data[,c(4:8)])
	
	#Normalize CA/RT values to fit on same plot
	CA_RT_norm = t(apply(CA_RT, 2, function(x) (x-min(x))/(max(x)-min(x))))+0.01
	
	#Define predictor names			   
	groups_CA_RT = factor(make.unique(colnames(data)[c(4:8)]), 
						 levels = make.unique(colnames(data)[c(4:8)]))
	
	#Create data track with normalized CA/RT values across window
	dtrack2 = DataTrack(data = CA_RT_norm, groups = groups_CA_RT, start = data$start,
						end = data$start+999, chromosome = chr, genome = gen, 
						name = "CA or RT/KB", type = c("smooth"), 
						span = 0.2)
	
	#Create needleplot with all tracks combined					 
	plotTracks(list(itrack, gtrack, grtrack, dtrack1, dtrack2), 
			   from =  data$start[1], to = data$start[nrow(data)],  
			   transcriptAnnotation = 'symbol', collapseTranscripts = 'meta', 
			   main = paste0(cohort_name, "\nObs=", obs, " Exp=", exp, " FDR=", FDR), 
			   background.title = "darkblue", cex.main = 0.9, extend.right = 20000, 
			   extend.left = 20000)
	
}

#Get paths for plot data					 
paths = list.files(pff(paste0("data/005H_needle_plot_data_ct/")), full.names = T)

#Create needle plots
pdf(pff("data/005H_needleplot_allct_2.pdf"), width  = 5, height = 6)
lapply(paths, plot_data)
dev.off()