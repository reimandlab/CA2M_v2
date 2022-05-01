source("000_HEADER.R")

source(pff("/bin/999_process_data.R"))

#Get data paths for ENCODE ATAC-Seq bigwigs
track_names = list.files(paste0(input_data_dir, "ENCODE_ATACSeq/"))[-c(1,2,208)]
track_paths = list.files(paste0(input_data_dir, "ENCODE_ATACSeq/"), 
						 full.names = T)[-c(1,2,208)]

#Get hg38 genome
genome = BSgenome.Hsapiens.UCSC.hg38

#Load in supplemetary file for the bigwigs
supp_file = fread(paste0(input_data_dir, 
						 "ENCODE_ATACSeq/experiment_report_2021_7_21_19h_35m.tsv"))[,c(2, 8)]

#Convert ATAC-seq bigwig to average per genomic window									 
get_ENCODE_ATACSeq_dataset = function(window_size){
	
    #Bin genome into windows
	gr.windows = tileGenome(seqinfo(Hsapiens)[paste0("chr", 1:22)], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom = TRUE)
	
    #Function to convert bigwig to binned average for each window
	merge_tracks = function(track_index){
		
		bigwig_path = list.files(paste0(track_paths[track_index], 
										"/released/GRCh38/fold_change_over_control/bigwig/rep1/"), 
								 full.names = T)[1] #Get bigwig path

		windowed_track = create_windowed_track(import(bigwig_path), window_size)[[3]] #Use function from script to convert to windowed track

	return(windowed_track)

    }
	
    #Merge windowed tracks from all samples into one data frame
    scores = as.data.table(do.call("cbind", lapply(seq(length(track_names)), 
													 merge_tracks)))
	#Add track names from supplementary file
    colnames(scores) = supp_file$`Biosample term name`[match(track_names, supp_file$Accession)]
	
	#Add genomic coordinates
    scores_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
												  start = start(gr.windows), 
												  scores))

    #Filter out windows with less than or equal to 80% mappability
    scores_dt_mappable = filter_windows_with_UMAP(scores_dt, window_size)
	
    #Make unique colnames as some are repeated
    colnames(scores_dt_mappable) = make.unique(colnames(scores_dt_mappable))

    return(scores_dt_mappable)}

#Create windowed tracks for different window sizes                                                                              
fwrite(get_ENCODE_ATACSeq_dataset(1000000), 
	   pff("/data/001D_ENCODE_ATACSeq_1MBwindow_processed.csv"))  
fwrite(get_ENCODE_ATACSeq_dataset(100000), 
	   pff("/data/001D_ENCODE_ATACSeq_100KBwindow_processed.csv"))                             
