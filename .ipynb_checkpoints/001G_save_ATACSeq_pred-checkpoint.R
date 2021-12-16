source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load in all ATAC-Seq and RT Preds
TCGA_ATAC = fread(pff("data/001C_TCGA_ATACSeq_1MBwindow_processed.csv"))
ENCODE_ATAC = fread(pff("data/001D_ENCODE_ATACSeq_1MBwindow_processed.csv"))
GEO_Brain_ATAC = fread(pff("data/001E_GEO_brain_ATACSeq_processed_1MB.csv"))
GEO_CLL_ATAC = fread(pff("data/001E_GEO_CLL_ATACSeq_processed_1MB.csv"))
GEO_HEK_ATAC = fread(pff("data/001E_GEO_HEK293_ATACSeq_processed_1MB.csv"))
GEO_lymphoma_ATAC = fread(pff("data/001E_GEO_Lymphoma_ATACSeq_processed_1MB.csv"))
GEO_NHM1_ATAC = fread(pff("data/001E_GEO_NHM1_ATACSeq_processed_1MB.csv"))
GEO_Prostate_ATAC = fread(pff("data/001E_GEO_Prostate_ATACSeq_processed_1MB.csv"))
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))

#Save all preds
All_preds_1MB = as.data.table(cbind.data.frame(chr = TCGA_ATAC$chr, start = TCGA_ATAC$start, TCGA_ATAC[,-c(1,2)],
											 ENCODE_ATAC, GEO_Brain_ATAC[,-c(1,2)], GEO_CLL_ATAC[,-c(1,2)], 
											  GEO_HEK_ATAC[,-c(1,2)], GEO_lymphoma_ATAC[,-c(1,2)], GEO_NHM1_ATAC[,-c(1,2)],
											 GEO_Prostate_ATAC[,-c(1,2)], RT[,-c(1,2)]))

dim(All_preds_1MB)
# [1] 2465  871

fwrite(All_preds_1MB, pff("data/001G_All_preds_1MB.csv"))


#Save tumor preds
PrimaryTumor_preds_1MB = as.data.table(cbind.data.frame(chr = TCGA_ATAC$chr, start = TCGA_ATAC$chr, TCGA_ATAC[,-c(1,2)],
														GEO_CLL_ATAC[,24:57], GEO_Prostate_ATAC[,7:12]))

dim(PrimaryTumor_preds_1MB)
# [1] 2465  423

fwrite(PrimaryTumor_preds_1MB, pff("data/001G_PrimaryTumor_preds_1MB.csv"))

#Save normal preds
Normal_preds_1MB = as.data.table(cbind.data.frame(chr = TCGA_ATAC$chr, start = TCGA_ATAC$chr, 
												  ENCODE_ATAC[,-c(1, 2, 5, 9, 14, 56, 92, 102, 109, 182, 183)], GEO_Brain_ATAC[,-c(1,2)], 
												  GEO_CLL_ATAC[,3:23], GEO_HEK_ATAC[,-c(1,2)], GEO_NHM1_ATAC[,-c(1,2)], 
												  GEO_Prostate_ATAC[,3:6]))

dim(Normal_preds_1MB)
# [1] 2465  343												  
fwrite(Normal_preds_1MB, pff("data/001G_Normal_preds_1MB.csv"))



#Load in all ATAC-Seq and RT Preds 100KB
TCGA_ATAC = fread(pff("data/001C_TCGA_ATACSeq_100KBwindow.csv"))
ENCODE_ATAC = fread(pff("data/001D_ENCODE_ATACSeq_100KBwindow_processed.csv"))
GEO_Brain_ATAC = fread(pff("data/001E_GEO_brain_ATACSeq_processed_100KB.csv"))
GEO_CLL_ATAC = fread(pff("data/001E_GEO_CLL_ATACSeq_processed_100KB.csv"))
GEO_HEK_ATAC = fread(pff("data/001E_GEO_HEK293_ATACSeq_processed_100KB.csv"))
GEO_lymphoma_ATAC = fread(pff("data/001E_GEO_Lymphoma_ATACSeq_processed_100KB.csv"))
GEO_NHM1_ATAC = fread(pff("data/001E_GEO_NHM1_ATACSeq_processed_100KB.csv"))
GEO_Prostate_ATAC = fread(pff("data/001E_GEO_Prostate_ATACSeq_processed_100KB.csv"))
RT = fread(pff("data/001F_ENCODE_repliseq_100KBwindow_processed.csv"))

#Save all preds
All_preds_100KB = as.data.table(cbind.data.frame(chr = TCGA_ATAC$chr, start = TCGA_ATAC$start, TCGA_ATAC[,-c(1,2)],
											 ENCODE_ATAC[,-c(1,2)], GEO_Brain_ATAC[,-c(1,2)], GEO_CLL_ATAC[,-c(1,2)], 
											  GEO_HEK_ATAC[,-c(1,2)], GEO_lymphoma_ATAC[,-c(1,2)], GEO_NHM1_ATAC[,-c(1,2)],
											 GEO_Prostate_ATAC[,-c(1,2)], RT[,-c(1,2)]))

dim(All_preds_100KB)
# [1] 24458   871

fwrite(All_preds_100KB, pff("data/001G_All_preds_100KB.csv"))