source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

source(pff("/bin/999_run_randomforest_experiment.R"))

#Load 1 MB mutation count dataset
mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_1MBwindow_processed.csv"))[ ,.SD, .SDcols = -c(1, 2)] 
                   
#Load in predictors
Normal_preds = fread(pff("data/001G_Normal_preds_1MB.csv"))
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]
Preds = as.data.table(cbind.data.frame(Normal_preds[,-c(1,2)], RT))

#Get observed vs predicted data tables
get_dt = function(cohort_index){
	
    cohort_name = colnames(mutation_counts_dt)[cohort_index]
    
    output = mutation_counts_dt[[cohort_index]]
    
    #Remove hypermutated windows from lymphoma and leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Normal_preds = Normal_preds[-c(292, 2438)]
        Preds = Preds[-c(292, 2438)]
        output = output[-c(292, 2438)]
    }
    
    rf = randomForest(x = Preds, 
					y = output, 
					keep.forest = T, 
					ntree = 1000, 
					do.trace = F, 
					importance = F)
    
    #Save observed vs. predicted values
    dt = as.data.table(cbind.data.frame(chr = Normal_preds$chr, 
										start = Normal_preds$start, 
										observed = output, 
										predicted = rf$predicted))
    
    fwrite(dt, paste0(pff("data/002E_normalCAplusRT_RF_obsvsexpected/"), cohort_name, ".csv"))}

mclapply(1:26, get_dt, mc.cores = 16)