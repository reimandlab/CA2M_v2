#Paths for scripts
home_dir = ""
setwd(home_dir)
input_data_dir = ""
PCAWG_MAF_path = ""
RepliSeq_path = ""
repliseq_supp_path = ""
Sig_annotation_file_path = "/.mounts/labs/reimandlab/private/generated_raw_data/PCAWG/PCAWG May2016 Data Release v1 %2F v1.1 %2F v1.2 %2F v1.3 %2F v1.4 - release_may2016v1.4 (1).csv"
PCAWG_Sigs_path = "/.mounts/labs/reimandlab/private/generated_raw_data/PCAWG/SigProfilier_PCAWG_WGS_probabilities_SBS.csv"

packages = c("data.table", "ActivePathways", "parallel", "gridExtra", "grid", 
			 "ggplot2", "mgsub", "gplots", "RColorBrewer", "GenomicRanges", 
			 "rtracklayer", "BSgenome.Hsapiens.UCSC.hg38", "plyr", "randomForest", 
			 "caret", "rCGH", "ggrastr", "cowplot", "pheatmap", "ggsci", "ggpubr",
			 "ggrepel", "gtools", "ggnewscale", "stringi" ,"Gviz")

lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE)) )

pff = function(x = "") {
	paste0("./", paste0(x, collapse = ""), sep = "")
}

opff = function(x) {
	system(paste("open", pff(x)))
}




file_open_call = function(fname) {
	fname1 = tail(strsplit(fname, '/')[[1]], 1)
	cat(paste0("\n\nscp -C $hn:", getwd(), "/", fname, " ./ && open ", fname1, "\n\n"))
}

file_open_call2 = function(fname) {
	system(paste("cp ", fname, "~/tmp/"))
	fname1 = tail(strsplit(fname, '/')[[1]], 1)
	cat(paste0("\n\n scp -C $cw:~/tmp/", fname1, " ./ && open ", fname1, "\n\n"))
}