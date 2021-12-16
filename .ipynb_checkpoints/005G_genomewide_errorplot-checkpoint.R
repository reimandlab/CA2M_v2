source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
library(Gviz)
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Load in qval error dt
qval_dt = fread(pff("data/005C_qvaldt_100KBerrorwindows.csv"))

#Define genome
gen = "hg38"

#Chromosomes
chrs = paste0("chr", 1:22)

#Function to plot errors across a chromosome
plot_errors = function(chr_in, cancer_type){
	 
	data = qval_dt[chr == chr_in]
	errors = -1*log10(data[[cancer_type]])
	
	itrack = IdeogramTrack(genome = gen, chromosome = chr_in)
	gtrack = GenomeAxisTrack()
	dtrack = DataTrack(data = errors, start = data$start,
						end = data$start+99999, chromosome = chr_in, genome = gen, 
						name = "-log10(FDR) of error", type = "p")
	
	plotTracks(list(itrack, gtrack, dtrack), 
			   from = data$start[1], 
			   to = data$start[nrow(data)]+99999)
	
}

pdf(pff("data/005F_needle_plots/melanoma_genomewide_errors.pdf"))
lapply(chrs, plot_errors, "Skin-Melanoma")
dev.off()


#Load in cytoband info
cyto_info = fread(paste0(input_data_dir, "cytoBandIdeo.txt"))
cyto_centromeres = cyto_info[V5 == "acen"]
cyto_centromeres.gr = GRanges(cyto_centromeres$V1, IRanges(cyto_centromeres$V2, cyto_centromeres$V3))

#See if high error windows are significantly enriched in centromere regions\
high_error_windows = qval_dt[`Skin-Melanoma` <= 10e-10]
high_error_windows.gr = GRanges(high_error_windows$chr, IRanges(high_error_windows$start, high_error_windows$start + 99999))
low_error_windows = qval_dt[`Skin-Melanoma` > 10e-10]
low_error_windows.gr = GRanges(low_error_windows$chr, IRanges(low_error_windows$start, low_error_windows$start + 99999))

#Get u-test comparing windows overlapping and not overlapping centromeres
higherror_centromere = length(unique(subjectHits(findOverlaps(cyto_centromeres.gr, high_error_windows.gr))))
higherror_nocentromere = nrow(high_error_windows) - higherror_centromere
lowerror_centromere = length(unique(subjectHits(findOverlaps(cyto_centromeres.gr, low_error_windows.gr))))
lowerror_nocentromere = nrow(low_error_windows) - lowerror_centromere

dat = data.frame(
  "high_error" = c(higherror_centromere, higherror_nocentromere),
  "low_error" = c(lowerror_centromere, lowerror_nocentromere),
  row.names = c("Overlapping_centromere", "Not_overlapping_centromere"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("high_error", "low_error")

test = fisher.test(dat)
expected_values = chisq.test(dat)$expected
#Only expects 0.6 windows to be high error and overlapping a centromeric region