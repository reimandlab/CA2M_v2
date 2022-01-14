source("000_HEADER.R")

#Get paths to RF model errors from 100KB-scale analysis
error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 
#Plot histograms of errors 
plot_hist = function(cohort_index){
    
	#Get project code
	project_code = error_dt_cancertypes[cohort_index]
    
	#Load in data table of observed and expected values
    error_dt = fread(error_dt_paths[cohort_index])
    
	#Define errors as observed - predicted mutations
    error_dt$errors = as.numeric(error_dt$observed - error_dt$predicted)
    
	#Create histogram of errors
    p = ggplot(error_dt, aes(x = errors)) + 
        geom_histogram(color = "black", fill = "steelblue", bins = 50)+
        labs(title = project_code, x = "Model error (observed-predicted mutations)", y = "Count")+
        theme_bw()
}

#Repeat for 23 cancer types									 
pdf(pff("data/005B_error_histograms.pdf"), width = 5, height = 5)
lapply(1:23, plot_hist)
dev.off()