source("000_HEADER.R")

source(pff("/bin/999_process_data.R"))

#Get table matching ATAC-Seq samples to TCGA samples
matching_table = fread(paste0(input_data_dir, "TCGA_ATACSEQ_identifier_mapping.txt"))

#Get unique TCGA Patients with ATAC seq
ATAC_patients = unique(unlist(lapply(matching_table$Case_ID, 
									 function(x) paste(unlist(strsplit(x, split = "-"))[1:3], collapse = "-"))))
ATAC_patients_all = unlist(lapply(matching_table$Case_ID, 
								 function(x) paste(unlist(strsplit(x, split = "-"))[1:3], collapse = "-")))
								  
#Load in TCGA tissue source site code table								  
TSS_codetable = fread(paste0(input_data_dir, "TCGA_tissueSourceSite.tsv"))

#Get abbreviations for each TCGA cancer type
TCGA_Abbreviations = fread(paste0(input_data_dir, "TCGA_diseaseStudy.tsv"))

#Get human genome hg38
genome = BSgenome.Hsapiens.UCSC.hg38

#Import ATAC-Seq bigwig paths and codes
bigwigpaths = list.files(paste0(input_data_dir,"TCGA_ATACSEQ_Individual_Bigwigs/"), full.names = T)
bigwig_patients_code = unlist(lapply(list.files(paste0(input_data_dir,"TCGA_ATACSEQ_Individual_Bigwigs/")),
	function(x) unlist(strsplit(x,split = ".", fixed = T))[1]))

#Convert ATAC-seq bigwig to average per genomic window									 
process_atacseq = function(window_size){                                 

    #Bin genome into windows
    gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr",1:22,sep="")], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom=TRUE)

    #Function to convert bigwig to binned average for each window
    merge_tracks = function(path){

        gr.data = import(path) #Load in bigwig

        windowed_track  =create_windowed_track(gr.data, window_size)[[3]] #Use function from script to convert to windowed track

        return(windowed_track)
    }

    #Merge windowed tracks from all samples into one data frame
    ATAC_SEQ_scores = do.call("cbind.data.frame", 
							  mclapply(bigwigpaths, 
									   merge_tracks, 
									   mc.cores = 16))
	#Add genomic coordinates
    ATAC_SEQ_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
												 start = start(gr.windows),
												 ATAC_SEQ_scores))
    
    #Process patient code names
    bigwig_patients_code = unlist(lapply(bigwig_patients_code, 
										 function(x) gsub("_","-",x)))
										 
    bigwig_patients_code_normal = ATAC_patients_all[match(bigwig_patients_code, matching_table$bam_prefix)]
    colnames(ATAC_SEQ_dt) = c("chr", "start", bigwig_patients_code_normal)

    ###Change column names to include cancer type as well as TCGA patient code		
										 
	#Get cancer types for each TCGA patient code									 
    ATACSeq_Sample_cancertypes = TSS_codetable$`Study Name`[match(unlist(lapply(colnames(ATAC_SEQ_dt)[-c(1, 2)], 
																				function(x) unlist(strsplit(x, split = "-", fixed = T))[2])), 
																				TSS_codetable$`TSS Code`)]
																		 
	#Convert TCGA cancer types to abbreviations																	 
    ATACSeq_Sample_abbreviations = TCGA_Abbreviations$`Study Abbreviation`[match(ATACSeq_Sample_cancertypes, 
																				 TCGA_Abbreviations$`Study Name`)]
																		 
	#Paste TCGA cancer type and TCGA patient code
    colnames(ATAC_SEQ_dt)[-c(1, 2)] = paste(ATACSeq_Sample_abbreviations, 
										 colnames(ATAC_SEQ_dt)[-c(1, 2)], 
										   sep = " ")

    #Filter out windows with less than or equal to 80% mappability
    ATAC_SEQ_dt_mappable = filter_windows_with_UMAP(ATAC_SEQ_dt, window_size)

    #Keep tracks from samples with more than 1 replicate
    cols_to_keep = names(table(colnames(ATAC_SEQ_dt)[-c(1,2)])[as.numeric(which(table(colnames(ATAC_SEQ_dt)[-c(1,2)])>1))])

    #Average tracks from same patient
    ATAC_SEQ_dt_mappable_norep=as.data.table( # sapply returns a list here, so we convert it to a data.frame
        sapply(cols_to_keep, # for each unique column name
           function(col) rowMeans(data.matrix(ATAC_SEQ_dt_mappable[,.SD,.SDcols=which(names(ATAC_SEQ_dt_mappable)==col)])) # calculate row means
        )
      )
	
	#Add genomic coordinates																	 
    ATAC_SEQ_dt_mappable_norep = as.data.table(cbind.data.frame(chr = ATAC_SEQ_dt_mappable$chr, 
																start = ATAC_SEQ_dt_mappable$start,
																ATAC_SEQ_dt_mappable_norep))
    
    #Remove LGG sample with positive correlation
    ATAC_SEQ_dt_final = ATAC_SEQ_dt_mappable_norep[,.SD,.SDcols = -which(colnames(ATAC_SEQ_dt_mappable_norep) == "LGG TCGA-FG-A4MY")]
                                                                              
    return(ATAC_SEQ_dt_final)}                                                                          
                                                                              

                                                                              
#Create windowed tracks for different window sizes                                                                              
fwrite(process_atacseq(1000000), pff("/data/001C_TCGA_ATACSeq_1MBwindow_processed.csv"))
fwrite(process_atacseq(100000), pff("/data/001C_TCGA_ATACSeq_100KBwindow.csv"))
