source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

#Load in sample size dataset
PCAWG_mutations_dt_hg38_SNP = fread(pff("data/001A_PCAWG_mutations_hg38.csv"))
PCAWG_sample_table = PCAWG_mutations_dt_hg38_SNP[, 
												 .(Cohort_size = uniqueN(Donor_ID)), 
												 by = Project_Code][order(Cohort_size, decreasing = T)]
n_mut_sig_project = PCAWG_mutations_dt_hg38_SNP[, .N, by = .(Project_Code)]

PCAWG_full_name = c("Liver Hepatocellular carcinoma", "Pancreas Adenocarcinoma", "Prostate adenocarcinoma",
				  "Breast Adenocarcinoma", "Kidney Renal cell carcinoma", "Medulloblastoma", 
				  "Ovary Adenocarcinoma",  "Diffuse large B-cell lymphoma", "Esophageal Adenocarcinoma", 
				  "Chronic lymphocytic leukemia", "Pilocytic Astrocytoma", "Pancreatic Neuroendocrine tumour", 
				  "Stomach Adenocarcinoma", "Melanoma", "Head/Neck Squamous cell carcinoma", 
				  "Thyroid Adenocarcinoma", "Lung Squamous cell carinoma",  "Kidney Renal cell carcinoma,\nchromophobe type",
				  "Colorectal Adenocarcinoma", "Uterus Adenocarcinoma", 
				  "Bone Osteosarcoma", "Glioblastoma multiforme","Bone Leiomyosarcoma", 
				  "Biliary Adenocarcinoma", "Lung Adenocarcinoma", "Bladder Transitional cell carcinoma", 
				  "Myeloid Myeloproliferative neoplasm", "Oligodendroglioma", "Cervix Squamous cell carcinoma",
                  "Myeloid Acute myeloid leukemia","Breast Lobular carcinoma", "Bone neoplasm, epithelioid", 
				  "Chondroblastoma", "Breast In situ adenocarcinoma", "Lymphoid (Not otherwise specified)", 
				  "Myelodysplastic syndrome", "Cervix Adenocarcinoma")

PCAWG_full = as.data.table(cbind.data.frame("Study Name" = PCAWG_full_name, 
											"Project Code" = PCAWG_sample_table$Project_Code, 
											"Number of samples" = PCAWG_sample_table$Cohort_size, 
											"Number of mutations" = n_mut_sig_project[match(PCAWG_sample_table$Project_Code, 
																							n_mut_sig_project$Project_Code)]$N)) 

PCAWG_full = as.data.table(rbind.data.frame(PCAWG_full, list("", "", "", ""), 
											list("TOTAL", "PANCAN", sum(PCAWG_full$`Number of samples`), 
												 sum(n_mut_sig_project$N))))

cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")

PCAWG_full$`Included/Excluded` = ifelse(PCAWG_full$`Project Code` %in% cancer_types_to_keep, "Included", "Excluded")     
PCAWG_full$`Included/Excluded`[c(38, 39)] = ""
                
pdf(pff("data/SF1_PCAWG_stat_table.pdf"), width = 12, height = 12)
grid.table(PCAWG_full, rows = NULL)
dev.off()