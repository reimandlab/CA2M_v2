#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

mutation_counts_dt = fread(pff("data/001B_PCAWG_mutation_counts_100KBwindow_processed.csv"))[,-c(1, 2)] 

#Load in predictors
Preds = fread(pff("data/001G_All_preds_100KB.csv"))

cohort_name = gsub(".csv", "", list.files(pff("data/001I_PCAWG_sigs_new/")))[-c(2,3,4,5,6,8,9,10,11,14,26,27)][cohort_index]

output = mutation_counts_dt[[cohort_name]]

#Remove hypermutated regions of the genome
Preds = Preds[-c(2858:2866, 20254, 24212:24219)]
output = output[-c(2858:2866, 20254, 24212:24219)]

dim(Preds[, -c(1,2)])
# [1] 24440   869

rf = randomForest(x = Preds[, -c(1,2)], 
				  y = output, 
				  keep.forest = T, 
				  ntree = 1000, 
				  do.trace = F, 
				  importance = F)

dt = as.data.table(cbind.data.frame(chr = Preds$chr,
									start = Preds$start, 
									observed = output, 
									predicted = rf$predicted))

fwrite(dt, paste0(pff("/data/005A_100KB_RF_errorDT_results/"), cohort_name, ".csv"))