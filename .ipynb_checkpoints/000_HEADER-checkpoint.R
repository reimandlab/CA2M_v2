home_dir = ""
setwd(home_dir)
input_data_dir = ""

PCAWG_MAF_path = ""

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