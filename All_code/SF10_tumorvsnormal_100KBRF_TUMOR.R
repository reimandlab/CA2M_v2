#Get cohort number as argument to loop through
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load megebase-scale binned mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 

#Core cancer types to keep
cancer_types_to_keep = c("PANCAN", "Breast-AdenoCa", "Prost-AdenoCA", 
						 "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", 
						 "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	

#Get cohort name
project_code = cancer_types_to_keep[cohort_index]

#Response of model will be mutation counts from cohort of interest
Output = mutation_counts_dt[[project_code]]                                    
 
#Load in all predictors at 100KB scale
Preds = fread(pff("data/001G_All_preds_100KB.csv"))

#Remove hypermutated windows if cohort is lymphoma or leukemia
if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
	Preds = Preds[-c(2858:2866, 20254, 24212:24219)]
	Output = Output[-c(2858:2866, 20254, 24212:24219)]
}

#Load in supplementary file for predictors
Pred_supp = fread(pff("data/001K_predictor_descriptions.csv"))           

#Combine primary cancer ATAC-Seq and replication timing predictors
Tumor_preds = Preds[,.SD,.SDcols = Pred_supp[Predictor_categories == "Primary cancer"]$Predictor_names]
RT = Preds[,.SD,.SDcols = Pred_supp[`Predictor_categories` == "RT"]$Predictor_names]
Tumor_RT_Preds = as.data.table(cbind.data.frame(Tumor_preds, RT))

#Run random forest with monte carlo cross-validation
RF_result_tumor = run_RF_and_get_pred_importances(Tumor_RT_Preds, Output, 100, n_tree = 1000, cores = 10, train_split = 0.8)

#Save results
fwrite(RF_result_tumor, paste0(pff("data/SF10_tumorvnormal_100KB/Tumor/"), project_code, ".csv"))

