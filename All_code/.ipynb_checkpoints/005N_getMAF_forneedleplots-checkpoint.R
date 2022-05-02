source("000_HEADER.R")

PCAWG_mutations_dt_hg38 = fread(pff("data/001A_PCAWG_mutations_hg38.csv"))

PCAWG_mutations_dt_hg38.gr = GRanges(PCAWG_mutations_dt_hg38$Chromosome, 
							  IRanges(PCAWG_mutations_dt_hg38$Start_position,
									  PCAWG_mutations_dt_hg38$End_position))

CGC_gene_windows.gr = readRDS(pff("data/005I_CGC_gene_window_gr.RDS"))

#Function to save MAF of specific gene and cancer type window
get_MAF = function(gene, cancer_type){
	gene_window.gr = CGC_gene_windows.gr[which(CGC_gene_windows.gr$gene == gene)]
	cancer_type_MAF = PCAWG_mutations_dt_hg38[Project_Code == cancer_type]
	cancer_type_MAF.gr = PCAWG_mutations_dt_hg38.gr[which(PCAWG_mutations_dt_hg38$Project_Code == cancer_type)]
	
	gene_cancer_type_window_MAF = cancer_type_MAF[subjectHits(findOverlaps(gene_window.gr, cancer_type_MAF.gr))]
	
	fwrite(gene_cancer_type_window_MAF, paste0(pff("data/005N_needleplot_MAFs/"), gene, "_", cancer_type, ".csv"))

}

get_MAF("PIK3CA", "Breast-AdenoCa")

get_MAF("IRF4", "Panc-AdenoCA")

get_MAF("HIF1A", "Thy-AdenoCA")