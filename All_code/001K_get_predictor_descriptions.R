source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load in all megabase-scale predictors
TCGA_ATAC = fread(pff("data/001C_TCGA_ATACSeq_1MBwindow_processed.csv"))[,-c(1,2)]
ENCODE_ATAC = fread(pff("data/001D_ENCODE_ATACSeq_1MBwindow_processed.csv"))
GEO_Brain_ATAC = fread(pff("data/001E_GEO_brain_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_CLL_ATAC = fread(pff("data/001E_GEO_CLL_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_HEK_ATAC = fread(pff("data/001E_GEO_HEK293_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_lymphoma_ATAC = fread(pff("data/001E_GEO_Lymphoma_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_NHM1_ATAC = fread(pff("data/001E_GEO_NHM1_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_Prostate_ATAC = fread(pff("data/001E_GEO_Prostate_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]

#Set descriptions and categories for each predictor in each dataset
TCGA_ATAC_descriptions = unlist(lapply(colnames(TCGA_ATAC), function(x) unlist(strsplit(x, split = " "))[1]))
TCGA_ATAC_categories = rep("Primary cancer", ncol(TCGA_ATAC))

ENCODE_ATAC_descriptions = colnames(ENCODE_ATAC)
ENCODE_ATAC_categories = rep("Normal tissue", ncol(ENCODE_ATAC))
ENCODE_ATAC_categories[c(5, 9, 14, 56, 92, 102, 109, 182, 183)] = "Cancer cell line"

GEO_Brain_descriptions = colnames(fread(paste0(input_data_dir, "Brain_GEO_descriptions.txt")))
GEO_Brain_categories = rep("Normal tissue", ncol(GEO_Brain_ATAC))

GEO_CLL_descriptions = unlist(lapply(strsplit(colnames(GEO_CLL_ATAC), split = "_"), 
								function(x) paste0(x[-c(1,2)], collapse = "_")))
GEO_CLL_descriptions = gsub(".bigWig", "", GEO_CLL_descriptions )
GEO_CLL_descriptions = gsub("_0d_", "_", GEO_CLL_descriptions )
GEO_CLL_categories = c(rep("Normal tissue", 21), rep("Primary cancer", 34))

GEO_HEK_descriptions = rep("HEK293", 4)
GEO_HEK_categories = rep("Normal tissue", 4)

GEO_lymphoma_descriptions = c("Raji", "Raji")
GEO_lymphoma_categories = c("Cancer cell line", "Cancer cell line")

GEO_NHM1_descriptions = "NHM1"
GEO_NHM1_categories = "Normal tissue"

GEO_Prostate_descriptions = c(rep("Prostate normal", 4), rep("Prostate cancer", 6))
GEO_Prostate_categories = c(rep("Normal tissue", 4), rep("Primary cancer", 6))

cell_phases = c("G1", "G2", "S1", "S2", "S3", "S4")
RT_Description = c(paste0("Embryonic stem cell ", cell_phases), 
				   paste0("Foreskin fibroblast ", rep(cell_phases, 2)), 
				   paste0("B-Lymphocyte ", rep(cell_phases, 5)), 
				   paste0("HeLa ", cell_phases), 
				   paste0("HepG2 ", cell_phases), 
				   paste0("HUVEC ", cell_phases), 
				   paste0("IMR90 ", cell_phases), 
				   paste0("K562 ", cell_phases), 
				   paste0("MCF-7 ", cell_phases), 
				   paste0("NHEK ", cell_phases), 
				   paste0("SK-N-SH ", cell_phases))
RT_categories = rep("RT", 96)

#Get names, descriptions, categories of all predictors
Predictor_names = c(colnames(TCGA_ATAC), colnames(ENCODE_ATAC), colnames(GEO_Brain_ATAC), colnames(GEO_CLL_ATAC),
					colnames(GEO_HEK_ATAC), colnames(GEO_lymphoma_ATAC), colnames(GEO_NHM1_ATAC),
					colnames(GEO_Prostate_ATAC), colnames(RT))
									 
Predictor_descriptions = c(TCGA_ATAC_descriptions, ENCODE_ATAC_descriptions, GEO_Brain_descriptions, GEO_CLL_descriptions,
							GEO_HEK_descriptions, GEO_lymphoma_descriptions, GEO_NHM1_descriptions, 
							GEO_Prostate_descriptions, RT_Description)

Predictor_categories = c(TCGA_ATAC_categories, ENCODE_ATAC_categories, GEO_Brain_categories, GEO_CLL_categories,
						GEO_HEK_categories, GEO_lymphoma_categories, GEO_NHM1_categories, GEO_Prostate_categories,
						RT_categories)

#Load dictionary of matching cancer type of each predictor
matching_predictor_dict = readRDS(pff("data/001J_Matching_predictors_dictionary.RDS"))		
									 
#Get matching cancer types for each predictor									 
Predictor_matching_PCAWG = unlist(lapply(Predictor_names, 
										 function(x) paste(names(matching_predictor_dict)[grep(x, matching_predictor_dict)], collapse = ", ")))

#Get original study where each predictor comes from										 
Predictor_study = c(rep("https://pubmed.ncbi.nlm.nih.gov/30361341/", ncol(TCGA_ATAC)), 
				    rep("https://pubmed.ncbi.nlm.nih.gov/32728249/", ncol(ENCODE_ATAC)), 
				    rep("https://pubmed.ncbi.nlm.nih.gov/29945882/", ncol(GEO_Brain_ATAC)), 
				    rep("https://pubmed.ncbi.nlm.nih.gov/31996669/", ncol(GEO_CLL_ATAC)), 
				    rep("https://pubmed.ncbi.nlm.nih.gov/30791920/", ncol(GEO_HEK_ATAC)), 
					rep("https://pubmed.ncbi.nlm.nih.gov/25753668/", ncol(GEO_lymphoma_ATAC)), 
					rep("https://pubmed.ncbi.nlm.nih.gov/29149598/", ncol(GEO_NHM1_ATAC)), 
					rep("https://pubmed.ncbi.nlm.nih.gov/32690948/", ncol(GEO_Prostate_ATAC)), 
					rep("https://www.pnas.org/content/107/1/139", ncol(RT)))

#Combine all information into one data table
predictor_dt = cbind.data.frame(Predictor_names, 
								Predictor_descriptions, 
								Predictor_categories, 
							    Predictor_matching_PCAWG, 
							    Predictor_study)
										 
#Save predictor supplementary data table										 
fwrite(predictor_dt, pff("data/001K_predictor_descriptions.csv"))