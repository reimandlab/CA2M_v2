source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

source(pff("/bin/999_process_data.R"))

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Import human genome build hg38                            
genome = BSgenome.Hsapiens.UCSC.hg38
            
#Import liftover files
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

#Function to get windowed track
process_ATACSeq_dataset = function(paths, names, window_size){
	
	gr.windows = tileGenome(seqinfo(Hsapiens)[paste0("chr", 1:22)], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom = TRUE)

	merge_tracks = function(path){
		print(which(paths == path))
		gr.data = import(path)

		#Liftover to hg38
		gr.data_hg38 = unlist(liftOver(gr.data, ch))
		
		windowed_track = create_windowed_track(gr.data_hg38, window_size)[[3]] #Use function from script to convert to windowed track

	return(windowed_track)

    }

    scores = as.data.table(do.call("cbind", mclapply(paths, merge_tracks, mc.cores = 4)))
    colnames(scores) = names
    scores_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
												  start = start(gr.windows), 
												  scores))

    #Filter out windows with less than or equal to 80% mappability
    scores_dt_mappable = filter_windows_with_UMAP(scores_dt, window_size)
	
    #Make unique colnames as some are repeated
    colnames(scores_dt_mappable) = make.unique(colnames(scores_dt_mappable))

    return(scores_dt_mappable)}


#Process brain tracks
brain_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Brain_tracks/"), full.names = T)
brain_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Brain_tracks/"))

brain_processed_1MB = process_ATACSeq_dataset(brain_paths, brain_names, 1000000)
fwrite(brain_processed_1MB, pff("data/001E_GEO_brain_ATACSeq_processed_1MB.csv"))

brain_processed_100KB = process_ATACSeq_dataset(brain_paths, brain_names, 100000)
fwrite(brain_processed_100KB, pff("data/001E_GEO_brain_ATACSeq_processed_100KB.csv"))

#Process CLL tracks
CLL_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/CLL_tracks/"), full.names = T)
CLL_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/CLL_tracks/"))

CLL_processed_1MB = process_ATACSeq_dataset(CLL_paths, CLL_names, 1000000)
fwrite(CLL_processed_1MB, pff("data/001E_GEO_CLL_ATACSeq_processed_1MB.csv"))

CLL_processed_100KB = process_ATACSeq_dataset(CLL_paths, CLL_names, 100000)
fwrite(CLL_processed_100KB, pff("data/001E_GEO_CLL_ATACSeq_processed_100KB.csv"))

#Process HEK293 tracks
HEK293_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/HEK293_tracks/"), full.names = T)
HEK293_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/HEK293_tracks/"))

HEK293_processed_1MB = process_ATACSeq_dataset(HEK293_paths, HEK293_names, 1000000)
fwrite(HEK293_processed_1MB, pff("data/001E_GEO_HEK293_ATACSeq_processed_1MB.csv"))

HEK293_processed_100KB = process_ATACSeq_dataset(HEK293_paths, HEK293_names, 100000)
fwrite(HEK293_processed_100KB, pff("data/001E_GEO_HEK293_ATACSeq_processed_100KB.csv"))

#Process lymphoma tracks
Lymphoma_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Lymphoma_tracks/"), full.names = T)
Lymphoma_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Lymphoma_tracks/"))

Lymphoma_processed_1MB = process_ATACSeq_dataset(Lymphoma_paths, Lymphoma_names, 1000000)
fwrite(Lymphoma_processed_1MB, pff("data/001E_GEO_Lymphoma_ATACSeq_processed_1MB.csv"))

Lymphoma_processed_100KB = process_ATACSeq_dataset(Lymphoma_paths, Lymphoma_names, 100000)
fwrite(Lymphoma_processed_100KB, pff("data/001E_GEO_Lymphoma_ATACSeq_processed_100KB.csv"))

#Process prostate tracks
Prostate_paths = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Prostate_tracks/"), full.names = T)
Prostate_names = list.files(paste0(input_data_dir, "GEO_ATACSeq_tracks/Prostate_tracks/"))

Prostate_processed_1MB = process_ATACSeq_dataset(Prostate_paths, Prostate_names, 1000000)
fwrite(Prostate_processed_1MB, pff("data/001E_GEO_Prostate_ATACSeq_processed_1MB.csv"))

Prostate_processed_100KB = process_ATACSeq_dataset(Prostate_paths, Prostate_names, 100000)
fwrite(Prostate_processed_100KB, pff("data/001E_GEO_Prostate_ATACSeq_processed_100KB.csv"))

#Process NHM1 tracks
NHM1_path = paste0(input_data_dir, "GEO_ATACSeq_tracks/GSM2476338_NHM1-ATAC.bw")

NHM1_processed_1MB = process_ATACSeq_dataset(NHM1_path, "GSM2476338_NHM1-ATAC.bw", 1000000)
fwrite(NHM1_processed_1MB, pff("data/001E_GEO_NHM1_ATACSeq_processed_1MB.csv"))

NHM1_processed_100KB = process_ATACSeq_dataset(NHM1_path, "GSM2476338_NHM1-ATAC.bw", 100000)
fwrite(NHM1_processed_100KB, pff("data/001E_GEO_NHM1_ATACSeq_processed_100KB.csv"))