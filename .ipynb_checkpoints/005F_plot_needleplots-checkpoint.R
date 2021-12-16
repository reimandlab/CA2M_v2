source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Load in gencode for gene track
GENCODE = fread(paste0(input_data_dir,"GENCODE_hg38_PROCESSED.txt"))[gene_type == "protein_coding"]
genes = TxDb.Hsapiens.UCSC.hg38.knownGene

plot_data = function(data){
	gen = "hg38"
	chr = data$chr[1]
	GENCODE = GENCODE[chr == chr]
	itrack = IdeogramTrack(genome = gen, chromosome = chr)
	gtrack = GenomeAxisTrack()
# 	grtrack = GeneRegionTrack(start = data$start[1], end = data$start[nrow(data)]+999,
# 							rstart=GENCODE$start, rends=GENCODE$end, chromosome=chr, genome=gen,
# 							symbol=GENCODE$gene_name,
# 							name="GENCODE", strand=GENCODE$strand)
	
# 	grtrack = BiomartGeneRegionTrack(start = data$start[1], end = data$start[nrow(data)]+999,
# 									 chromosome = chr, genome = gen)
# 	grtrack <- GeneRegionTrack(genes, genome = gen,
#                            chromosome = chr, name = "UCSC known genes",
#                            exonAnnotation='symbol')
	grtrack = GeneRegionTrack(paste0(input_data_dir, "gencode.v38.annotation.gtf"),
							  stacking = 'squish', stackHeight = 0.3, name = 'Genes')
	dtrack1 = DataTrack(data = data$mut, start = data$start,
						end = data$start+999, chromosome = chr, genome = gen, 
						name = "Mutations/KB", type = c("p", "smooth"))
	dtrack2 = DataTrack(data = data$ATAC, start = data$start,
						end = data$start+999, chromosome = chr, genome = gen, 
						name = "CA/KB", type = c("p", "smooth"))
	
	plotTracks(list(itrack, gtrack, grtrack, dtrack1, dtrack2), 
			   from =  data$start[1]-20000, to = data$start[nrow(data)]+999+20000,   
			   transcriptAnnotation = 'symbol', collapseTranscripts = 'longest')
	
}

pdf(pff("data/005F_needle_plots/14_cancer_type_errorwindow.pdf"))
plot_data(fread(pff("data/005E_needle_plot_data/14_cancer_type_window.csv")))
dev.off()

pdf(pff("data/005F_needle_plots/13_cancer_type_errorwindow.pdf"))
plot_data(fread(pff("data/005E_needle_plot_data/13_cancer_type_window.csv")))
dev.off()

pdf(pff("data/005F_needle_plots/10_cancer_type_errorwindow_2.pdf"))
plot_data(fread(pff("data/005E_needle_plot_data/10_cancer_type_window.csv")))
dev.off()