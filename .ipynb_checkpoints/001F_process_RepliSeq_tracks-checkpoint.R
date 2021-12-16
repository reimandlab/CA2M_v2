date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("/bin/999_process_data.R"))
	   
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

RepliSeq_paths = list.files("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2MutData/RepliSeq_Tracks/",full.names=T)[-c(97,98)]
RepliSeq_names = unlist(lapply(list.files("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2MutData/RepliSeq_Tracks/"),
	function(x) unlist(strsplit(x,split=".bigWig"))[1]))[-c(97,98)]

#Import human genome build hg38                            
genome = BSgenome.Hsapiens.UCSC.hg38
                            
#Import liftover files
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

#Function to average individual tracks                           
merge_tracks = function(bigwig_path, window_size){

	gr.data = import(bigwig_path)

	#Liftover to hg38
	gr.data_hg38 = unlist(liftOver(gr.data, ch))

    windowed_track = create_windowed_track(gr.data_hg38, window_size)[[3]] #Use function from script to convert to windowed track
	
    return(windowed_track)

}                            
   

#Function to process and merge tracks into one data table                            
process_repliseq_data = function(paths, names, window_size){
    
    gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr",1:22,sep = "")], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom = TRUE)
    
    scores= as.data.table(do.call("cbind",mclapply(paths, 
													function(x) merge_tracks(x, window_size), 
													mc.cores = 8)))
    colnames(scores) = names
    scores_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
											   start = start(gr.windows), 
											   scores))
    
    scores_mappable = filter_windows_with_UMAP(scores_dt, window_size)
	
	#Get supplementary data
	supp = fread("/.mounts/labs/reimandlab/private/users/oocsenas/ChrAcc2MutData/repliseq_supplementary.tsv")
	supp_file_names = unlist(lapply(supp$Files, function(x) unlist(strsplit(x,split="/"))[3]))

	colnames(scores_mappable)[-c(1,2)] = make.names(unlist(lapply(colnames(scores_mappable)[-c(1,2)], 
										function(x) supp[grep(x, supp$Files)]$Description)), 
										unique = T)

# 	colnames(scores_mappable)[-c(1,2)] = gsub("Repli.seq.of.", "", colnames(scores_mappable)[-c(1,2)])                                                  
	colnames(scores_mappable)[-c(1,2)] = gsub("\\.", " ", colnames(scores_mappable)[-c(1,2)])  
    
    return(scores_mappable)
}      
                                                                               
fwrite(process_repliseq_data(RepliSeq_paths, RepliSeq_names, 100000), 
	   pff("/data/001F_ENCODE_repliseq_100KBwindow_processed.csv"))                             
fwrite(process_repliseq_data(RepliSeq_paths, RepliSeq_names, 1000000), 
	   pff("/data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))                             
                             
                             
                             
                                                  