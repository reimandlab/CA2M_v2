source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

library(chemometrics)
error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                           
                                   
get_pvals = function(cohort_index, log){
    print(cohort_index)
    project_code = error_dt_cancertypes[cohort_index]
    
    error_dt = fread(error_dt_paths[cohort_index])
    
    error_dt$errors = as.numeric(error_dt$observed - error_dt$predicted)
	#Trim extreme values for mean and sd
	p_values = pnorm(error_dt$errors, 
					    mean = mean(error_dt$errors),
						sd = sd(error_dt$errors),
						lower.tail = F, 
						log.p = log) #Natural log

    if(cohort_index == 1){
        return(cbind.data.frame(chr = error_dt$chr, 
								start = error_dt$start, 
								p_values))
    }else{
        return(p_values)
    }
}

#Save p-val of errors									 
p_val_dt = as.data.table(do.call("cbind.data.frame", lapply(c(1:23), get_pvals, log = F)))                                   
colnames(p_val_dt)[-c(1, 2)] = error_dt_cancertypes[c(1:23)] 
fwrite(p_val_dt, pff("data/005C_pvaldt_100KBerrorwindows.csv"))									 								 
									 
#Save log p-val of errors									 
log_p_val_dt = as.data.table(do.call("cbind.data.frame", lapply(c(1:23), get_pvals, log = T)))                                   
colnames(log_p_val_dt)[-c(1, 2)] = error_dt_cancertypes[c(1:23)]                                   
fwrite(log_p_val_dt, pff("data/005C_logpvaldt_100KBerrorwindows.csv"))

#Save q-val of errors
q_val_dt = as.data.table(cbind.data.frame(p_val_dt[,c(1, 2)], 
										  apply(p_val_dt[,-c(1, 2)], 
												2, 
												p.adjust, 
												method = "fdr")))         									 
fwrite(q_val_dt, pff("data/005C_qvaldt_100KBerrorwindows.csv"))

#Get dt for 17 core cancer types
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
									 
qval_dt_core = q_val_dt[,.SD,.SDcols = cancer_types_to_keep]

#Get number of unique windows significant in at least one cancer type
sum(rowSums(qval_dt_core < 0.05, na.rm = T) > 0, na.rm = T)
# 1570

#Get number of windows significant for each cancer type (overlapping)
num_sig = apply(qval_dt_core, 2, function(x) sum(x < 0.05, na.rm = T))
num_sig				
#   Breast-AdenoCa    Prost-AdenoCA       Kidney-RCC    Skin-Melanoma 
#               61              121               28               93 
#   Uterus-AdenoCA      Eso-AdenoCa  Stomach-AdenoCA          CNS-GBM 
#               77              151              133              221 
#         Lung-SCC ColoRect-AdenoCA  Biliary-AdenoCA         Head-SCC 
#              167              123              109              110 
#        Lymph-CLL     Lung-AdenoCA       Lymph-BNHL        Liver-HCC 
#              109              152               53               39 
#      Thy-AdenoCA 
#              103 
				
#For each window get number of significant cancer types
num_sig_cancertypes = table(apply(qval_dt_core, 1, function(x) sum(x < 0.05, na.rm = T)))
								  
num_sig_cancertypes								  
#     0     1     2     3     4    12    14 
# 22870  1370   158    30     9     2     1 