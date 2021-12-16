source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Import liftover files
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

#Load in paths to ATAC_Seq bigwigs
ENCODE_ATAC_track_names = colnames(fread(pff("data/001D_ENCODE_ATACSeq_100KBwindow_processed.csv"))[,-c(1,2)])
ENCODE_ATAC_track_paths = list.files(paste0(input_data_dir, "ENCODE_ATACSeq/"), 
						 full.names = T)[-c(1,2,208)]

TCGA_ATAC_track_paths = list.files("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/TCGA_ATACSEQ_Individual_Bigwigs/", 
								   full.names=T)
TCGA_ATAC_track_names = unlist(lapply(list.files("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/TCGA_ATACSEQ_Individual_Bigwigs/"),
	function(x) unlist(strsplit(x,split = ".", fixed = T))[1]))

#Load in MAF file
PCAWG_mutations_dt_hg38 = fread(pff("data/001A_PCAWG_mutations_hg38.csv"))
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))[,-c(1, 2)][-c(2858:2866, 20254, 24212:24219)]

#Load in window qvals
q_val_dt = fread(pff("data/005C_qvaldt_100KBerrorwindows.csv"))

#Load in human genome hg38
genome = BSgenome.Hsapiens.UCSC.hg38

#Function to get mutations every 100 bp for 100kb window
get_100bp_data = function(window_index, cohort_name, predictor_track, build){
	
	#Get 100bp windows inside 100kb window 
	gr.window = GRanges(q_val_dt$chr[window_index], 
						IRanges(q_val_dt$start[window_index], q_val_dt$start[window_index] + 99999))

    gr.windows_small = tileGenome(seqinfo(Hsapiens)[q_val_dt$chr[window_index]], 
							tilewidth = 1000, 
							cut.last.tile.in.chrom=TRUE)
	
	gr.windows = gr.windows_small[subjectHits(findOverlaps(gr.window,gr.windows_small))]
	
	#Get mutations for each 100bp window
	
	#Load in mutation data
    mut = PCAWG_mutations_dt_hg38[Project_Code == cohort_name]
    
    #Convert mutations file into GRanges object
    mutation.gr = GRanges(mut$Chromosome, 
						  IRanges(mut$Start_position,mut$End_position))
    
    #Get number of counts of mutations overlapping each genomic window
    mutation_counts = as.numeric(table(factor(queryHits(findOverlaps(gr.windows,mutation.gr)), 
											  levels = 1:length(gr.windows))))
	
	#Get ATAC-seq for each 100bp window
	
	#Import bigwig
	BigWig = import(predictor_track)
	
	#Liftover to hg38 if data is hg19
	if(build == "hg19"){
		BigWig = unlist(liftOver(BigWig, ch))
	}
	
	#Get an RleList object with one coverage vector per seqlevel in the BigWig
    gr.data.new = GenomicRanges::coverage(BigWig, weight="score")
    

	
    #Match seqlevels between input and human genome
    seqlevels(gr.windows, pruning.mode="coarse") = names(gr.data.new)
    
    #Get binned average for each window
    gr.data.new = binnedAverage(gr.windows, gr.data.new, "value")
	
	
	final_dt = cbind.data.frame(chr = as.character(seqnames(gr.windows)),
								start = start(gr.windows),
								mut = mutation_counts,
							    ATAC = gr.data.new$value)
    
	return(final_dt)
}


#Get tracks for highest error windows
# liver_error = which.min(q_val_dt$`Liver-HCC`)
liver_error = 8123
liver_track = list.files(paste0(ENCODE_ATAC_track_paths[which(ENCODE_ATAC_track_names == "left lobe of liver")], 
				  "/released/GRCh38/fold_change_over_control/bigwig/rep1/"), full.names = T)[1]
window_data = get_100bp_data(liver_error, "Liver-HCC", liver_track)
fwrite(window_data, pff("data/005E_needle_plot_data/liver_top.csv"))

#Plot first most common error gene (14 cancer types)
# skin_error = which.min(q_val_dt$`Skin-Melanoma`)
top_window = 17620								  
top_track = TCGA_ATAC_track_paths[662]
window_data = get_100bp_data(top_window, "Skin-Melanoma", top_track, "hg38")
fwrite(window_data, pff("data/005E_needle_plot_data/14_cancer_type_window.csv"))

#Plot second and third most commonly high error window
top_window = 23093								  
top_track = TCGA_ATAC_track_paths[662]
window_data = get_100bp_data(top_window, "Skin-Melanoma", top_track, "hg38")
fwrite(window_data, pff("data/005E_needle_plot_data/13_cancer_type_window.csv"))

									  
top_window = 9555								  
top_track = TCGA_ATAC_track_paths[662]
window_data = get_100bp_data(top_window, "Skin-Melanoma", top_track, "hg38")
fwrite(window_data, pff("data/005E_needle_plot_data/10_cancer_type_window.csv"))
