source("000_HEADER.R")

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 
                   
#Load in primary cancer CA and RT predictors
Tumor_preds = fread(pff("data/001G_PrimaryTumor_preds_1MB.csv"))
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]
Preds = as.data.table(cbind.data.frame(Tumor_preds[,-c(1,2)], RT))

#Get observed vs predicted data tables
get_dt = function(cohort_index){
	
	#Get name of cohort of interest
    cohort_name = colnames(mutation_counts_dt)[cohort_index]
    
	#Output of RF model is binned mutations counts from cohort of interest
    output = mutation_counts_dt[[cohort_index]]
    
    #Remove hypermutated windows if cohort is lymphoma or leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Tumor_preds = Tumor_preds[-c(292,2438)]
        Preds = Preds[-c(292,2438)]
        output = output[-c(292,2438)]
    }
    
	#Train random forest to predict binned mutation counts using primary cancer CA + RT
    rf = randomForest(x = Preds, 
					y = output, 
					keep.forest = T, 
					ntree = 1000, 
					do.trace = F, 
					importance = F)
    
    #Save observed vs. predicted values
    dt = as.data.table(cbind.data.frame(chr = Tumor_preds$chr, 
										start = Tumor_preds$start, 
										observed = output, 
										predicted = rf$predicted))
    
    fwrite(dt, paste0(pff("data/002D_tumourCAplusRT_RF_obsvsexpected/"), cohort_name, ".csv"))}

#Apply through all PCAWG cohorts
mclapply(1:26, get_dt, mc.cores = 16)