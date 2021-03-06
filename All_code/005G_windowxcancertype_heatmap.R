source("000_HEADER.R")

#Load in data table with q-values for 100KB windows
window_qval_dt = fread(pff("data/005C_qvaldt_100KBerrorwindows.csv"))

#Get corresponding genomic ranges for 100KB windows
windows_gr = GRanges(window_qval_dt$chr, 
					 IRanges(window_qval_dt$start, 
							 window_qval_dt$start+99999))

#Only keep Chromosomes 1-22,X in GENCODE
chr_to_keep = paste("chr", 
					c(1:22,"X"), 
					sep = "")

#Load in gene data table from GENCODE
GENCODE = fread(paste0(input_data_dir, "GENCODE_hg38_PROCESSED.txt"))[gene_type == "protein_coding"][chr %in% chr_to_keep]

#Convert genomic coordinates of genes to genomic ranges
GENCODE_gr = GRanges(GENCODE$chr, IRanges(GENCODE$start, GENCODE$end))

#Keep only core 17 cancer types for window FDR table						   
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")		
window_qval_dt_core = window_qval_dt[,.SD,.SDcols = cancer_types_to_keep]

#Keep windows which are significant in at least 1/17 cancer types
window_qval_dt_core_sig = window_qval_dt_core[which(rowSums(window_qval_dt_core < 0.05, na.rm = T) > 0)]
window_sig_gr = windows_gr[which(rowSums(window_qval_dt_core < 0.05, na.rm = T) > 0)]

#Load in known cancer genes from the CGC dataset
CGC = fread(paste0(input_data_dir, "Census_allFri Mar 26 16 35 45 2021.csv"))
CGC_genes = CGC$`Gene Symbol`

#Get most overlapping significant window for each gene
get_window = function(gene_index){
	
	#Get gene name
	gene_name = GENCODE[gene_index]$gene_name
	
	#Get gene genomic coordnates for gene
	gene_coordinates = GENCODE_gr[gene_index]

	#Get overlapping windows for gene
	overlapping_windows = subjectHits(findOverlaps(gene_coordinates, window_sig_gr))
	
	if(length(overlapping_windows) == 0){ #If no overlapping significant windows then return NA
            window_index = NA
        }

    if(length(overlapping_windows) == 1){   #If only 1 overlapping significant windows then return that one
            window_index = overlapping_windows
        }
        
    if(length(overlapping_windows) > 1){ #If multiple overlapping significant windows then return one with most overlap
            
            #Get amount of overlap for each overlapping significant window
            widths = unlist(lapply(overlapping_windows, 
								   function(x) width(intersect(gene_coordinates, 
															   window_sig_gr[x]))))
            
			#Get window with most overlap
            top_window = overlapping_windows[which.max(widths)]
            
            window_index = top_window                   
        }
        
    return(window_index)
}

#This list contains index of overlapping window for each gene
gene_windows = unlist(mclapply(1:nrow(GENCODE), get_window, mc.cores = 16))
								 
#Get overlapping genes for each window of interest
get_overlapping_genes = function(window_num){
	
	#Get genomic ranges of window
    window.gr = window_sig_gr[window_num]

    #Get overlapping genes for window using previous list
    overlapping_genes = which(gene_windows == window_num)
    overlapping_gene_names = GENCODE$gene_name[overlapping_genes]
    
    #Keep only known cancer genes from CGC
    overlapping_gene_names = overlapping_gene_names[which(overlapping_gene_names %in% CGC_genes)]
    
    if(length(overlapping_gene_names) >0){
        #Comma seperate all overlapping genes
        gene_list=paste(overlapping_gene_names,collapse=", ")}
	
    else{
        gene_list = ""
    }

return(gene_list)}

#List contains overlapping CGC genes for windows of interest
gene_list = unlist(lapply(1:nrow(window_qval_dt_core_sig), get_overlapping_genes))

#Remove 0 p-values and convert to data matrix for FDR of windows of interest
window_qval_dt_core_sig_nozero = as.data.table(apply(window_qval_dt_core_sig, 2, 
												   function(x) ifelse(x < 10**-20, 10**-20, x)))
qval_matrix = -1*log10(t(data.matrix(window_qval_dt_core_sig_nozero)))

#Create color palette
my_palette = colorRampPalette(c("white", "red"))(n = 299)

#Plot heatmap of significant window FDR values for the different cancer types
#Windows are annotated by overlapping CGC genes
pdf(pff("/data/005I_windowxcancertype_heatmap2_cordist.pdf"), width = 50, height = 5)
heatmap.2(qval_matrix, margins = c(12,8), col = my_palette, 
		  distfun = function(x) as.dist(sqrt(2*(1-cor(t(x))))), 
		  cexRow = 0.8, cexCol = 0.2, lhei = c(1,4), 
          lwid = c(1, 20), keysize = 2, key.par = list(cex = 0.5), 
		  labCol = gene_list, density.info = "none", 
          trace = "none")
dev.off()
 
#Get which windows overlap CGC genes in each cancer type and the genes that overlap them
CGC_gene_windows = gene_windows[which(GENCODE$gene_name %in% CGC_genes)]
CGC_gene_windows.gr = window_sig_gr[na.omit(CGC_gene_windows)]
CGC_gene_windows.gr$gene = gene_list[na.omit(CGC_gene_windows)]
saveRDS(CGC_gene_windows.gr, pff("data/005I_CGC_gene_window_gr.RDS"))				 

###Are significant windows more likely to overlap known cancer genes than other genes?		

#Get all genes overlapping any window													 
all_genes_mappable = GENCODE$gene_name[unique(queryHits(findOverlaps(GENCODE_gr, windows_gr)))]

#Get genes which overlap a significantly high error window
genes_overlapping_higherror = GENCODE$gene_name[unique(queryHits(findOverlaps(GENCODE_gr, window_sig_gr)))]
length(genes_overlapping_higherror)													 
# [1] 900

#Get genes which are known cancer genes and overlap any window											 
CGC_genes_mappable = CGC_genes[which(CGC_genes %in% all_genes_mappable)]

#Exact fisher's test find enrichment of known cancer genes among genes overlapping high error windows
p_val = fisher.test(all_genes_mappable %in% CGC_genes_mappable, 
					all_genes_mappable %in% genes_overlapping_higherror,
					alternative = "g")$p.val
p_val												 
# [1] 3.108516e-08		

#Expected number of CGC genes overlapping high error windows
Expected_value = (length(CGC_genes_mappable) * length(genes_overlapping_higherror))/ length(all_genes_mappable)
Expected_value																				  
#[1] 33.21731
													 
#Observed number of CGC genes overlapping high error windows
CGC_overlapping_genes = intersect(CGC_genes_mappable, genes_overlapping_higherror)
length(CGC_overlapping_genes)
# [1] 67

###Create table of known cancer genes overlapping high error windows and their corresponding information
CGC_overlapping_genes_cts = unlist(lapply(CGC_overlapping_genes, 
										 function(x) paste(colnames(window_qval_dt_core_sig)[which(window_qval_dt_core_sig[gene_windows[which(GENCODE$gene_name == x)]] < 0.05)], collapse = ", ")))													 
CGC_gene_CGC_cts = CGC$`Tumour Types(Somatic)`[match(CGC_overlapping_genes, CGC$`Gene Symbol`)]
										  
table_dt = as.data.table(cbind.data.frame(Gene = CGC_overlapping_genes, 
										  Cancer_types = CGC_overlapping_genes_cts, 
										  CGC_cancer_types = CGC_gene_CGC_cts))										  

table_dt$Cancer_types = unlist(lapply(table_dt$Cancer_types, 
									 function(x) paste(strwrap(x, width = 30), collapse = "\n")))
table_dt$CGC_cancer_types = unlist(lapply(table_dt$CGC_cancer_types, 
									 function(x) paste(strwrap(x, width = 30), collapse = "\n")))		
										  
table_dt$Role = CGC$`Role in Cancer`[match(CGC_overlapping_genes, CGC$`Gene Symbol`)]										  
										  
library(gridExtra)
library(grid)										  
table = tableGrob(table_dt, rows = NULL, 
					theme = ttheme_default(base_size = 6))	
									 									 
									 
									 
fname = pff("data/005I_significant_CGCgene.pdf")			  
pdf(fname, 
		width = 10, height = 20)
grid.draw(table)
dev.off()	
										  
										  
