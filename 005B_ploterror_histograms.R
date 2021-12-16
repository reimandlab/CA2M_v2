source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 

plot_hist = function(cohort_index){
    print(cohort_index)
    project_code = error_dt_cancertypes[cohort_index]
    
    error_dt = fread(error_dt_paths[cohort_index])
    
    error_dt$errors = as.numeric(error_dt$observed - error_dt$predicted)
    
    p = ggplot(error_dt, aes(x = errors)) + 
        geom_histogram(color = "black", fill = "steelblue", bins = 50)+
        labs(title = project_code, x = "Model error (observed-predicted mutations)", y = "Count")+
        theme_bw()
}
                                   
pdf(pff("data/005B_error_histograms.pdf"),width=5,height=5)
lapply(1:23, plot_hist)
dev.off()