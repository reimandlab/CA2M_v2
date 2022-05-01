source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load in cohort names in analysis
cohort_names = unlist(lapply(list.files(pff("/data/002G_permuted_normalvtumor_results/Tumor_results/")), 
										function(x) unlist(strsplit(x, 
																	split = ".RDS", 
																	fixed = T))[1]))

#Load in tumor CA + all RT random forest results
data_paths_tumor_RT = list.files(pff("/data/002G_permuted_normalvtumor_results/Tumor_results/"), full.names = T)
accuracies_tumor_RT = unlist(lapply(data_paths_tumor_RT, 
								 function(x) readRDS(x)))


#Load in normal CA + all RT random forest results
data_paths_normal_RT = list.files(pff("/data/002G_permuted_normalvtumor_results/Normal_results/"), full.names = T)
accuracies_normal_RT = unlist(lapply(data_paths_normal_RT, 
								 function(x) readRDS(x)))

#Combine results into one data table
accuracies.dt = as.data.table(cbind.data.frame(cohort_name = unlist(lapply(cohort_names, rep, 1000)), 
											   Tumor_adjR2 = accuracies_tumor_RT, 
											   Normal_adjR2 = accuracies_normal_RT))
									 
#Get the difference in accuracy between model trained on primary cancer CA + RT vs. normal tissue CA + RT								  
accuracies.dt$diff_adjR2 = accuracies.dt$Tumor_adjR2 - accuracies.dt$Normal_adjR2
                               
#Get median accuracies and confidence intervals into one table
plot_dt = accuracies.dt[,.(median_tumor_adjR2 = median(Tumor_adjR2), 
						   median_normal_adjR2 = median(Normal_adjR2), 
						   mean_diff_adjR2 = mean(diff_adjR2), 
						   low_CI_adjR2 = quantile(diff_adjR2, 0.025), 
						   high_CI_adjR2 = quantile(diff_adjR2, 0.975)), 
						by = cohort_name]                               
                
order = plot_dt$cohort_name[order(plot_dt$mean_diff_adjR2)]

#Get empirical p-value for each cancer type
get_emp_pval = function(cohort_name_in){
	data = accuracies.dt[cohort_name == cohort_name_in]
	mean_diff = mean(data$diff_adjR2)
	p_val  = sum(sign(data$diff_adjR2) != sign(mean_diff))/1000
	return(p_val)
	
}
emp_pvals = sapply(plot_dt$cohort_name, get_emp_pval)			   
			   
plot_dt$emp_pvals = emp_pvals
plot_dt$emp_FDR = p.adjust(emp_pvals, method="fdr")							  

plot_dt$bar_label = ifelse(plot_dt$emp_pvals < 0.001, "***",
						ifelse(plot_dt$emp_pvals < 0.01, "**",
							ifelse(plot_dt$emp_pvals < 0.05, "*",
								ifelse(plot_dt$emp_pvals < 0.1, "\u2022", ""))))		
plot_dt$score_position = unlist(lapply( plot_dt$low_CI_adjR2, function(x) min(0, x) - 0.02))

#Create barplot
cairo_pdf(pff("data/002H_CA_tumourvsnormal_permuted_barplots.pdf"), width = 9, height = 5)                       
ggplot(plot_dt, 
	aes(x = factor(cohort_name, 
					levels = rev(order)), 
					y = mean_diff_adjR2, 
					fill = median_tumor_adjR2)) + 
	geom_bar(stat = "identity", colour = 'black')+ #Add barplot
	geom_errorbar(aes(ymin = low_CI_adjR2, ymax = high_CI_adjR2), width = .2)+ #Add errorbars
	geom_text(aes(label = bar_label, y = high_CI_adjR2 + 0.02), fontface = "bold", size = 3)+ #Add significance symbol
	geom_text(aes(label = round(median_tumor_adjR2, 2), #Add text showing tumor model accuracy
				y = score_position), 
			size = 2, fontface = "bold")+	
	theme_bw()+
    scale_fill_gradientn(colours = c("white", "red"), limits = c(0, 1))+ #Add barplot fill 
	labs(y = "Change in model accuracy,\ntumour CA+RTvs. normal tissue CA+RT (adjusted R2)", #Add labels
		 x = "Cohort name", 
		 fill = "Model accuracy using\ntumour CA+RT\n(adjusted R2)")+
    theme(axis.title = element_text(size=9), #Format text
        axis.text.x = element_text(size = 8, 
								   angle = 90, 
								   hjust = 1, 
								   vjust = 0.5,
								   colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
         legend.text = element_text(size = 8),
         legend.title =element_text(size = 10))							   
dev.off()                               
							  
							  
#Create and save CSV's of statistical results
pvals = sapply(unique(accuracies.dt$cohort_name), 
			   function(x) wilcox.test(accuracies.dt[cohort_name==x]$diff_adjR2)$p.val)	
plot_dt$p_val = as.numeric(pvals[match(plot_dt$cohort_name, names(pvals))])
plot_dt$FDR = p.adjust(plot_dt$p_val, method="fdr")
	   
fwrite(accuracies.dt, pff("data/002H_resultsdata.csv"))
fwrite(plot_dt, pff("data/002H_results_summary.csv"))			   
