source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

mut_cor = cor(mutation_counts_dt, method = "spearman")
pdf(pff("/data/SF4_1MBmuttrack_cormap.pdf"), width = 10, height = 10)
heatmap.2(x = mut_cor, 
    Colv = T,
	Rowv = T,	  
    dendrogram = "both",
    col = "bluered",
    trace = "none",
    xlab = "Cancer Type",
	ylab = "Cancer Type",
	margins = c(13,13),
	cexRow = 1.3, cexCol = 1.3)
dev.off()