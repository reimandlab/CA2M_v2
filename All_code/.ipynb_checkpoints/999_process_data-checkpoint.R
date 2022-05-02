
create_windowed_track = function(BigWig, window_size){
    
    #'Take in a BigWig file for a whole-genome ChIP-Seq/DNase-Seq/ATAC-Seq/etc. track and convert it into genome-wide windows 
    #'of specified size with an average read count for each window.
    
    #'BigWig object imported using the import function. Should be GRanges with an additional score column.
    
    #'window_size can be any integer between 2 and 50,000,000
    
    #Load in human genome hg38
    genome = BSgenome.Hsapiens.UCSC.hg38
    
    #Convert genome (chromosomes 1 through 22) into windows of specified size and cut last window in each chromosome to stay within chromosome
    gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr",1:22,sep="")], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom=TRUE)
    
    #Get an RleList object with one coverage vector per seqlevel in the BigWig
    gr.data.new = GenomicRanges::coverage(BigWig, weight="score")
    
    #Match seqlevels between input and human genome
    seqlevels(gr.windows, pruning.mode="coarse") = names(gr.data.new)
    
    #Get binned average for each window
    gr.data.new = binnedAverage(gr.windows, gr.data.new, "value")
    
    #Convert to data table with chromosome and start base pairs of window in addition to the binned average value
    track_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)),
											  start = start(gr.windows), 
											  mean_accessibility = gr.data.new$value))

    return(track_dt)
}


create_RMV = function(MAF, window_size){
        
    #'Convert MAF mutations file (data.table or dataframe) into a regional mutational variation (RMV) track. 
    #'This means sum of mutations for each window throughout the genome for the specified window size
    
    #'MAF file can be a data.table or a dataframe with any columns. Required columns include columns with names of Chromosome, Start_position, and End_position.
    
    
    #Load in human genome hg38
    genome = BSgenome.Hsapiens.UCSC.hg38
    
    #Convert genome (chromosomes 1 through 22) into windows of specified size and cut last window in each chromosome to stay within chromosome
    gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr", 1:22, sep = "")], 
							tilewidth = window_size, 
							cut.last.tile.in.chrom=TRUE)
    
    #Convert mutations file into GRanges object
    mutation.gr = GRanges(MAF$Chromosome, 
						  IRanges(MAF$Start_position,MAF$End_position))
    
    #Get number of counts of mutations overlapping each genomic window
    mutation_counts = as.numeric(table(factor(queryHits(findOverlaps(gr.windows,mutation.gr)), 
											  levels = 1:length(gr.windows))))
    
    #Convert to data table with chromosome and start base pairs of window in addition to the mutation counts
    RMV_dt = as.data.table(cbind.data.frame(chr = as.character(seqnames(gr.windows)), 
											start = start(gr.windows), 
											Mut_counts = mutation_counts))
 
    return(RMV_dt)
}



filter_windows_with_UMAP=function(track, window_size){
    
    #' Takes in any windowed track (binnedAverage or RMV) and removes any windows which don't meet a mappability cutoff based on UMAP.
    #' Cutoff is that 80% of the window must be mappable.
    #' Takes about 75s right now. 
    
    #Convert track to GenomicRanges object
    track.gr = GRanges(track$chr, 
					   IRanges(track$start, track$start+window_size-1))
    
    #Load in and process UMAP mappability file
    UMAP = suppressWarnings(fread(""))
    colnames(UMAP)[1:3] = c("chr","start","end")
	
    #Get regions defined mappable and give them score of 1 to use with binned average
    UMAP.gr = disjoin(GRanges(UMAP$chr, IRanges(UMAP$start,UMAP$end))) #takes long
    UMAP.gr = GRanges(UMAP.gr, value = rep(1, length(UMAP.gr)))
    
    #Get an RleList object with one coverage vector per seqlevel in UMAP
    gr.umap.new = GenomicRanges::coverage(UMAP.gr, weight="value")
    
    #Match seqlevels between UMAP and input track
    seqlevels(track.gr, pruning.mode="coarse") = names(gr.umap.new)
    
    #Get binned average of UMAP file (value indicates fractiion of base pairs which are mappable for each window)
    gr.umap.new = binnedAverage(track.gr, gr.umap.new, "value")
    
    #Keep windows which have 80% or greater mappable base pairs 
    track_mappable = track[which(gr.umap.new$value > 0.8)]
    
    return(track_mappable)
    }

