source("000_HEADER.R")

source(pff("/bin/999_process_data.R"))
	   
#Load in Repli-seq track paths and names
RepliSeq_paths = list.files(RepliSeq_path, full.names = T)[-c(97,98)]
RepliSeq_names = unlist(lapply(list.files(RepliSeq_path),
	function(x) unlist(strsplit(x,split=".bigWig"))[1]))[-c(97,98)]

#Import human genome build hg38                            
genome = BSgenome.Hsapiens.UCSC.hg38
                            
#Get supplementary data
supp = fread(repliseq_supp_path)
supp_file_names = unlist(lapply(supp$Files, function(x) unlist(strsplit(x,split="/"))[3]))

#Function to average individual tracks                           
merge_tracks = function(bigwig_path, window_size){
	
	#Import bigwig
	gr.data = import(bigwig_path)

	#Liftover to hg38
	gr.data_hg38 = unlist(liftOver(gr.data, ch))

    windowed_track = create_windowed_track(gr.data_hg38, window_size)[[3]] #Use function from script to convert to windowed track
	
    return(windowed_track)

}                            
   

#Function to process and merge tracks into one data table                            
process_repliseq_data = function(paths, names, window_size){
    
	#Bin genome into windows
    gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr",1:22,sep = "")], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom = TRUE)
    
	#Merge windowed tracks from all samples into one data frame
    scores= as.data.table(do.call("cbind",mclapply(paths, 
													function(x) merge_tracks(x, window_size), 
													mc.cores = 8)))
    colnames(scores) = names
	
	#Add genomic coordinates
    scores_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
											   start = start(gr.windows), 
											   scores))
    #Filter out windows with less than or equal to 80% mappability
    scores_mappable = filter_windows_with_UMAP(scores_dt, window_size)
	

	#Add track name descriptions to column names
	colnames(scores_mappable)[-c(1,2)] = make.names(unlist(lapply(colnames(scores_mappable)[-c(1,2)], 
										function(x) supp[grep(x, supp$Files)]$Description)), 
										unique = T)

# 	colnames(scores_mappable)[-c(1,2)] = gsub("Repli.seq.of.", "", colnames(scores_mappable)[-c(1,2)])                                                  
	colnames(scores_mappable)[-c(1,2)] = gsub("\\.", " ", colnames(scores_mappable)[-c(1,2)])  
    
    return(scores_mappable)
}  
														   
#Create windowed tracks for different window sizes                                                                              
fwrite(process_repliseq_data(RepliSeq_paths, RepliSeq_names, 100000), 
	   pff("/data/001F_ENCODE_repliseq_100KBwindow_processed.csv"))                             
fwrite(process_repliseq_data(RepliSeq_paths, RepliSeq_names, 1000000), 
	   pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))                             
                             
                             
                             
                                                  