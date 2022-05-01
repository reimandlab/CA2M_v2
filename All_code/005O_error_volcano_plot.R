source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Load in qval dt
window_qval_dt = fread(pff("data/005C_qvaldt_100KBerrorwindows.csv"))

error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 

get_errors = function(cohort_name){
    print(cohort_name)
	cohort_index = which(error_dt_cancertypes == cohort_name)
	    
    error_dt = fread(error_dt_paths[cohort_index])
    
    error_dt$errors = as.numeric(error_dt$observed - error_dt$predicted)
	
	return(error_dt$errors)}
									 
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", 
						 "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")
#CANCER TYPES DONT MATCH UP									 
error_dt = as.data.table(do.call("cbind.data.frame", 
								 lapply(cancer_types_to_keep, get_errors)))			
colnames(error_dt) = cancer_types_to_keep
error_dt = as.data.table(cbind.data.frame(chr = window_qval_dt$chr,
													  start = window_qval_dt$start,
													  error_dt))									 
error.m = melt(error_dt, id.vars = c("chr", "start"))									 

window_qval_dt_ct = window_qval_dt[, .SD, .SDcols = c("chr", "start", cancer_types_to_keep)]
									 
window_qval_dt.m = melt(window_qval_dt_ct, id.vars = c("chr", "start"))
									 
final_dt = merge(error.m, window_qval_dt.m, by = c("chr", "start", "variable") )
colnames(final_dt)[c(3, 4, 5)] = c("cancer_type", "error", "qval")
									 
#Add labels of cancer genes
CGC = fread(paste0(input_data_dir, "Census_allFri Mar 26 16 35 45 2021.csv"))
CGC_genes = CGC$`Gene Symbol`
									 
GENCODE = fread(paste0(input_data_dir, "GENCODE_hg38_PROCESSED.txt"))[gene_type == "protein_coding"]
GENCODE_CGC = GENCODE[gene_name %in% CGC_genes]									 
GENCODE_CGC_gr = GRanges(GENCODE_CGC$chr, IRanges(GENCODE_CGC$start, GENCODE_CGC$end))
									 
final_dt.gr = GRanges(final_dt$chr, IRanges(final_dt$start, final_dt$start + 99999))

overlaps = findOverlaps(final_dt.gr, GENCODE_CGC_gr)								 
genes = unlist(lapply(unique(queryHits(overlaps)), 
			   function(x) paste(GENCODE_CGC$gene_name[subjectHits(overlaps)[which(queryHits(overlaps) == x)]], collapse = ", ")))	
	
final_dt$CGC_gene = rep("", nrow(final_dt))
final_dt$CGC_gene[unique(queryHits(overlaps))] = genes					  
					  
#Remove non-significant points
final_dt_sig = final_dt[qval<0.05]
					  
#Cap at q = 10e-50 and error = 500				  
final_dt_sig$qval[which(final_dt_sig$qval < 1e-50)] = 1e-50
final_dt_sig$error[which(final_dt_sig$error > 400)] = 400
					  
#Keep label if point is more significant than 1e-10
final_dt_sig$CGC_gene[which(final_dt_sig$qval > 1e-10 & !(final_dt_sig$CGC_gene %in% c("IRF4", "PIK3CA") & final_dt_sig$cancer_type  == "Breast-AdenoCa"))] = ""	
					  
#Add color column
PCAWG_colours = readRDS(paste0(input_data_dir, "PCAWG_colour_palette.RDS"))
# levels(final_dt_sig$cancer_type) = sort(levels(final_dt_sig$cancer_type))
				
final_dt_sig$cancer_type = as.character(final_dt_sig$cancer_type)				  

colors = as.character(PCAWG_colours)[match(tolower(sort(unique(final_dt_sig$cancer_type))), 
											  names(PCAWG_colours))]
final_dt_sig$cancer_type = factor(final_dt_sig$cancer_type, 
								  levels = sort(unique(final_dt_sig$cancer_type)))			  
					  
#Create volcano plot
pdf(pff("data/005O_error_volcanoplot.pdf"), height = 10, width = 12)
ggplot(final_dt_sig, aes(x = error, 
					 y = -1*log10(qval), 
					 fill = cancer_type,
					label = CGC_gene))+
	geom_point(pch = 21, size = 3, color = "black")+
	ylim(c(0, 50))+
	xlim(c(0, 400))+
	theme_bw()+
	theme(axis.text = element_text(size = 10))+
	scale_fill_manual(values = colors, labels = levels(final_dt_sig$cancer_type))+
	labs(x = "Residuals (observed mutations - predicted mutations)", 
		y = "-log10(FDR)", fill = "Cancer Type")+
	geom_text_repel(size = 4)
									 
dev.off()