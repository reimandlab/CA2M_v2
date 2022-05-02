source("000_HEADER.R")

error_dt_paths = list.files(pff("/data/005A_100KB_RF_errorDT_results"), full.names = T)
error_dt_cancertypes = unlist(lapply(list.files(pff("/data/005A_100KB_RF_errorDT_results")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))									 
									 
#Get mutations in IRF4 window for each cancer type
get_muts = function(cohort_name){
	
    
    error_dt = fread(error_dt_paths[which(error_dt_cancertypes == cohort_name)])
	
	window_index = which(error_dt$chr == "chr6" & error_dt$start == 300001)
	
	
	muts = error_dt$observed[seq(window_index-10, window_index+10)]
	predicted = error_dt$predicted[seq(window_index-10, window_index+10)]
	label = c(rep("background", 10), "interest", rep("background", 10))
	
	dt = cbind.data.frame(window_num = 1:21, 
						  muts,
						  predicted,
						  label, 
						  cohort_name = rep(cohort_name, 21))
	
	return(dt)
	
}									 
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")										 
plot_dt = as.data.table(do.call("rbind.data.frame", lapply(cancer_types_to_keep, get_muts)))	
		
									 
p = ggplot(plot_dt, aes(x = window_num,
							   y = muts))+
	theme_bw()+
	geom_bar(stat = "identity", aes(fill = label))+ 
	geom_line(stat = "identity", aes(y = predicted), color = "blue")+
	geom_point(stat = "identity", aes(y = predicted), color = "blue", size = 0.5)+								 
	facet_wrap(~cohort_name, ncol = 3, scales = "free")+
	scale_fill_manual(aes(y = muts), values = c("black", "red"), labels = c("background", "chr6:300,001:400,000"))+
	scale_color_manual(aes(y = predicted), name = "predicted")

									 
									 
pdf(pff("data/SF9_IRF4_barplots.pdf"), width = 8, height = 10)
p
dev.off()
									 
									 