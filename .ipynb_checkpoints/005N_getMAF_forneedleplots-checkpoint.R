source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

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
	
	fwrite(gene_cancer_type_window_MAF, paste0(pff("data/005N_needleplot_MAFs/", gene, "_", cancer_type, ".csv")))

}



#1 PIK3CA Breast
PIK3CA_window.gr = CGC_gene_windows.gr[which(CGC_gene_windows.gr$gene == "PIK3CA")]
Breast_MAF = PCAWG_mutations_dt_hg38[Project_Code == "Breast-AdenoCa"]
Breast_MAF.gr = PCAWG_mutations_dt_hg38.gr[which(PCAWG_mutations_dt_hg38$Project_Code == "Breast-AdenoCa")]

PIK3CA_breast_window_MAF = Breast_MAF[subjectHits(findOverlaps(PIK3CA_window.gr, Breast_MAF.gr))]
fwrite(PIK3CA_breast_window_MAF, pff("data/005N_needleplot_MAFs/PIK3CA_Breast_MAF.csv"))

#1 IRF4 Panc
IRF4_window.gr = CGC_gene_windows.gr[which(CGC_gene_windows.gr$gene == "IRF4")]
Breast_MAF = PCAWG_mutations_dt_hg38[Project_Code == "Breast-AdenoCa"]
Breast_MAF.gr = PCAWG_mutations_dt_hg38.gr[which(PCAWG_mutations_dt_hg38$Project_Code == "Breast-AdenoCa")]

PIK3CA_breast_window_MAF = Breast_MAF[subjectHits(findOverlaps(PIK3CA_window.gr, Breast_MAF.gr))]
fwrite(PIK3CA_breast_window_MAF, pff("data/005N_needleplot_MAFs/PIK3CA_Breast_MAF.csv"))