source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")
source(pff("/bin/999_run_RF_withSHAP.R"))

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load 1 MB binned mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

#Get SHAP values in each window for each predictor in each cancer type
get_SHAP_values = function(cohort_name){

	#Output for RF is binned mutation count in cohort
    output = mutation_counts_dt[[cohort_name]]

    #Remove hypermutated windows if cohort is lymphoma or leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
        output = output[-c(292, 2438)]
    }
    

    #Train randomForest and get SHAP values
    SHAP_values = run_RF_and_get_SHAP(Preds, output)
    
    fwrite(SHAP_values, paste0(pff("data/003F_RF_SHAP_Results/"), cohort_name, ".csv"))
}
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	


#Get SHAP values in all cohrot and save results
mclapply(cancer_types_to_keep, get_SHAP_values, mc.cores = 9)



                                                                 
