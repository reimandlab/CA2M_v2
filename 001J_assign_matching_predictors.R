source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load in predictors
TCGA_ATAC = fread(pff("data/001C_TCGA_ATACSeq_1MBwindow_processed.csv"))[,-c(1,2)]
ENCODE_ATAC = fread(pff("data/001D_ENCODE_ATACSeq_1MBwindow_processed.csv"))[,-c(1,2)]
GEO_Brain_ATAC = fread(pff("data/001E_GEO_brain_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_CLL_ATAC = fread(pff("data/001E_GEO_CLL_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_HEK_ATAC = fread(pff("data/001E_GEO_HEK293_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_lymphoma_ATAC = fread(pff("data/001E_GEO_Lymphoma_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_NHM1_ATAC = fread(pff("data/001E_GEO_NHM1_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
GEO_Prostate_ATAC = fread(pff("data/001E_GEO_Prostate_ATACSeq_processed_1MB.csv"))[,-c(1,2)]
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]


`PANCAN` = ""           

`Biliary-AdenoCA` = c(colnames(TCGA_ATAC)[221:237],
					colnames(ENCODE_ATAC)[grep("liver", colnames(ENCODE_ATAC))],
					 colnames(RT)[55:60])

`Bone-Leiomyo` = ""  

`Bone-Osteosarc` = ""

`Breast-AdenoCa` = c(colnames(TCGA_ATAC)[c(20:85)],
					 colnames(ENCODE_ATAC)[c(90, 144, 173, 197)],
					 colnames(RT)[79:84])

`CNS-GBM` = c(colnames(TCGA_ATAC)[c(145:153, 211:220)], 
			colnames(GEO_Brain_ATAC))      

`CNS-Medullo` = colnames(GEO_Brain_ATAC) 

`CNS-PiloAstro` = c(colnames(TCGA_ATAC)[c(145:153, 211:220)], 
			colnames(GEO_Brain_ATAC))    

`ColoRect-AdenoCA` = c(colnames(TCGA_ATAC)[93:129],
					colnames(ENCODE_ATAC)[grep("colon", colnames(ENCODE_ATAC))])

`Eso-AdenoCa` = c(colnames(TCGA_ATAC)[130:144], colnames(ENCODE_ATAC)[c(51, 56)])

`Head-SCC` = c(colnames(TCGA_ATAC)[c(154:160)],
				colnames(ENCODE_ATAC)[56],
				colnames(RT)[85:90]) 

`Kidney-ChRCC` = c(colnames(TCGA_ATAC)[c(161:210)], colnames(GEO_HEK_ATAC))     

`Kidney-RCC`= c(colnames(TCGA_ATAC)[c(161:210)], colnames(GEO_HEK_ATAC)) 

`Liver-HCC` = c(colnames(TCGA_ATAC)[c(221:237)],
				colnames(ENCODE_ATAC)[grep("liver", colnames(ENCODE_ATAC))], 
				colnames(RT)[55:60]) 

`Lung-AdenoCA` = c(colnames(TCGA_ATAC)[c(238:259)],
					colnames(ENCODE_ATAC)[grep("lung", colnames(ENCODE_ATAC))]) 

`Lung-SCC` = c(colnames(TCGA_ATAC)[c(260:275)],
					colnames(ENCODE_ATAC)[grep("lung", colnames(ENCODE_ATAC))],
				colnames(RT)[85:90])     

`Lymph-BNHL` = c(colnames(ENCODE_ATAC)[grep("K562", colnames(ENCODE_ATAC))],
				 colnames(ENCODE_ATAC)[grep("B cell", colnames(ENCODE_ATAC))],
				 colnames(GEO_lymphoma_ATAC), colnames(GEO_CLL_ATAC), colnames(RT)[19:48])

`Lymph-CLL` =  c(colnames(ENCODE_ATAC)[grep("K562", colnames(ENCODE_ATAC))],
				 colnames(ENCODE_ATAC)[grep("B cell", colnames(ENCODE_ATAC))],
				 colnames(GEO_lymphoma_ATAC), colnames(GEO_CLL_ATAC), colnames(RT)[19:48])

`Ovary-AdenoCA` = colnames(ENCODE_ATAC)[c(42,39, 79,203)]

`Panc-AdenoCA` = colnames(ENCODE_ATAC)[grep("pancreas", colnames(ENCODE_ATAC))] 

`Panc-Endocrine` = colnames(ENCODE_ATAC)[grep("pancreas", colnames(ENCODE_ATAC))]

`Prost-AdenoCA` = c(colnames(TCGA_ATAC)[c(291:316)], colnames(GEO_Prostate_ATAC)) 

`Skin-Melanoma` = c(colnames(TCGA_ATAC)[c(317:327)], colnames(GEO_NHM1_ATAC)) 

`Stomach-AdenoCA`= c(colnames(TCGA_ATAC)[c(328:347)], colnames(ENCODE_ATAC)[c(51, 56, 66, 146, 174)])

`Thy-AdenoCA` = c(colnames(TCGA_ATAC)[c(357:370)], colnames(ENCODE_ATAC)[grep("thyroid", colnames(ENCODE_ATAC))])  

`Uterus-AdenoCA` = c(colnames(TCGA_ATAC)[c(371:381)], colnames(ENCODE_ATAC)[c(39, 162)])  

Matching_preds = list(`PANCAN`=`PANCAN`,`Biliary-AdenoCA`=`Biliary-AdenoCA`,`Bone-Leiomyo`=`Bone-Leiomyo`,`Bone-Osteosarc`=`Bone-Osteosarc`,`Breast-AdenoCa`=`Breast-AdenoCa`, `CNS-GBM`=`CNS-GBM`,`CNS-Medullo`=`CNS-Medullo`,`CNS-PiloAstro`=`CNS-PiloAstro`,`ColoRect-AdenoCA`=`ColoRect-AdenoCA`,`Eso-AdenoCa`=`Eso-AdenoCa`, `Head-SCC`=`Head-SCC`,`Kidney-ChRCC`=`Kidney-ChRCC`,`Kidney-RCC`=`Kidney-RCC`,`Liver-HCC`=`Liver-HCC`,`Lung-AdenoCA`=`Lung-AdenoCA`,`Lung-SCC`=`Lung-SCC`,`Lymph-BNHL`=`Lymph-BNHL`,`Lymph-CLL`=`Lymph-CLL`,`Ovary-AdenoCA`=`Ovary-AdenoCA`,`Panc-AdenoCA`=`Panc-AdenoCA`,`Panc-Endocrine`=`Panc-Endocrine`,`Prost-AdenoCA`=`Prost-AdenoCA`,`Skin-Melanoma`=`Skin-Melanoma`,`Stomach-AdenoCA`=`Stomach-AdenoCA`,`Thy-AdenoCA`=`Thy-AdenoCA`,`Uterus-AdenoCA`=`Uterus-AdenoCA`)

saveRDS(Matching_preds, pff("data/001J_Matching_predictors_dictionary.RDS"))
