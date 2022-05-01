#Get cohort number and block as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])
block = as.integer(args[2])

source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

#Load 1 MB binned mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[,-c(1, 2)]

#Get name of cohort based on index
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
project_code = cancer_types_to_keep[cohort_index] 

#Output of RF model is binned mutatation count of cohort
output = mutation_counts_dt[[project_code]]

#Remove hypermutated windows if cohort is lymphoma or leukemia
if(project_code %in% c("Lymph-CLL", "Lymph-BNHL")){
    Preds = Preds[-c(292, 2438)]
    output = output[-c(292, 2438)]
}

#Permute model output to run permutation test based on random seed
get_RF_dt = function(seed){
    print(seed)
    set.seed(seed)

    #Randomly permute response without replacement
    output_scrambled = sample(as.numeric(output))

    #Train randomForest to predict scrambled output
    rf=randomForest(x = Preds, y = output_scrambled, keep.forest = T, 
					ntree = 1000, do.trace = F, importance = T)
    
    return(as.numeric(importance(rf, type = 1, scale = F)))}

#Blocks of seeds to iterate through to get 1000 permutations in total
blocks = seq(100, 1000, 100)           

#Run code on desired block using input argument
RF_results = as.data.table(do.call("rbind.data.frame", 
								   mclapply(seq(blocks[block]-99, blocks[block]), 
											get_RF_dt, mc.cores = 10)))
colnames(RF_results) = colnames(Preds)

#Save results of permutation test
dir.create(paste0(pff("/data/003B_Feature_importance_permutations/"), project_code, "/"))

fwrite(RF_results, paste0(pff("/data/003B_Feature_importance_permutations/"), project_code, "/", block, ".csv"))

print("Finished")
                                                                