source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

filepaths = list.files(pff("data/004A_Sig_RF_Results"), 
					   full.names = T)
cohort_names = unlist(lapply(list.files(pff("data/004A_Sig_RF_Results")),
                            function(x) unlist(strsplit(x ,split = ".csv"))[1]))			 
		 
get_plot_data = function(cohort_name){
    print(cohort_name)

    data = fread(filepaths[which(cohort_names == cohort_name)])
                            
    Sigs = fread(paste0(pff("data/001I_PCAWG_sigs_new/"), cohort_name, ".csv"))[,-c(1,2)]
    Sigs_to_keep =  colnames(Sigs[,.SD, .SDcols = which(colSums(Sigs)>20000)])

    Adj_R2 = as.numeric(data[1,.SD,.SDcols = Sigs_to_keep])

    get_mut_num = function(SBS){
		n_mut = sum(Sigs[[SBS]])
		return(n_mut)}
    mut_nums = unlist(lapply(Sigs_to_keep, get_mut_num))
	
    plot_dt = cbind.data.frame(cohort_name = rep(cohort_name, length(Sigs_to_keep)), 
							   signature = Sigs_to_keep, 
							   mut_nums, 
							   Adj_R2)
    
    return(plot_dt)
                            
                            
}

cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")

full_plot_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(cancer_types_to_keep, get_plot_data)))  
                            		 
etiology_dt = as.data.table(cbind.data.frame(
    etiology = c(rep("APOBEC enzyme activity", 2), 
				 rep("Defective DNA repair", 4), 
				 rep("SBS1", 1), 
				 rep("Exogenous/Carcinogen", 8),
				 rep("SBS5, SBS40", 2),
				 rep("Unknown/Other", 8)), 
    signature = c(c("SBS13", "SBS2"), 
				  c("SBS3", "SBS6", "SBS26", "SBS44"), 
				  c("SBS1"), 
				  c("SBS4", "SBS7a", "SBS7b", "SBS7c","SBS7d","SBS29", "SBS22", "SBS35"),
				  c("SBS5", "SBS40"),
				  c( "SBS9", "SBS12","SBS16", "SBS17a", "SBS17b", "SBS18", "SBS37", "SBS41"))))							 
                            
full_plot_dt$etiologies = etiology_dt$etiology[match(full_plot_dt$signature, 
											  etiology_dt$signature)]
                            
#Create boxplots
full_plot_dt$etiologies = factor(full_plot_dt$etiologies, 
								 levels = c("SBS1", 
											"APOBEC enzyme activity", 
											"Defective DNA repair", 
											"Exogenous/Carcinogen",
											"SBS5, SBS40",
											"Unknown/Other"))                   
           

endogenous_etiologies = c("APOBEC enzyme activity", "Defective DNA repair", "SBS1")							 
p_val_1 = wilcox.test(full_plot_dt[etiologies %in% endogenous_etiologies]$Adj_R2, 
					  full_plot_dt[etiologies == "Exogenous/Carcinogen"]$Adj_R2)$p.val
p_val_2 = wilcox.test(full_plot_dt[etiologies %in% endogenous_etiologies]$Adj_R2, 
					  full_plot_dt[etiologies == "SBS5, SBS40"]$Adj_R2)$p.val	
p_val_3 = wilcox.test(full_plot_dt[etiologies %in% endogenous_etiologies]$Adj_R2, 
					  full_plot_dt[etiologies == "Unknown/Other"]$Adj_R2)$p.val	
							 
							 
							 
#Run ANOVA tests to account for mutation burden of signatures
full_plot_dt$class = factor(ifelse(full_plot_dt$etiologies %in% endogenous_etiologies, 
								   "endogenous", as.character(full_plot_dt$etiologies)))

#Run for carcinogenic vs. endogenous signatures							 
H0 = lm(Adj_R2 ~ mut_nums, 
		data = full_plot_dt[class %in% c("endogenous", "Exogenous/Carcinogen")])
							 
H1 = lm(Adj_R2 ~ mut_nums + class, 
		data = full_plot_dt[class %in% c("endogenous", "Exogenous/Carcinogen")])	
							 
pval_anova_1 = anova(H0, H1)$`Pr(>F)`[2]							 
							 							 
#Run for SBS5/40 vs. endogenous signatures							 
H0 = lm(Adj_R2 ~ mut_nums, 
		data = full_plot_dt[class %in% c("endogenous", "SBS5, SBS40")])
							 
H1 = lm(Adj_R2 ~ mut_nums + class, 
		data = full_plot_dt[class %in% c("endogenous", "SBS5, SBS40")])	
							 
pval_anova_2 = anova(H0, H1)$`Pr(>F)`[2]							 
							 
#Run for unknown vs. endogenous signatures							 
H0 = lm(Adj_R2 ~ mut_nums, 
		data = full_plot_dt[class %in% c("endogenous", "Unknown/Other")])
							 
H1 = lm(Adj_R2 ~ mut_nums + class, 
		data = full_plot_dt[class %in% c("endogenous", "Unknown/Other")])	
							 
pval_anova_3 = anova(H0, H1)$`Pr(>F)`[2]									 
							 
#Keep certain labels
labels_to_keep = c("Skin-Melanoma\nSBS7a", 
				   "Lung-AdenoCA\nSBS4", 
				   "Lung-AdenoCA\nSBS40", 
				   "Breast-AdenoCa\nSBS5", 
				   "Breast-AdenoCa\nSBS1", 
				   "Breast-AdenoCa\nSBS13", 
				   "Breast-AdenoCa\nSBS3", 
				   "Eso-AdenoCa\nSBS17b", 
				   "CNS-GBM\nSBS40", 
				   "CNS-GBM\nSBS1", 
				   "Kidney-RCC\nSBS1", 
				   "Kidney-RCC\nSBS5", 
				   "Kidney-RCC\nSBS22",
				   "ColoRect-AdenoCA\nSBS40",
				   "Lung-SCC\nSBS13",
				   "Eso-AdenoCa\nSBS1",
				   "Stomach-AdenoCA\nSBS3", 
				   "Liver-HCC\nSBS4",
				   "Lung-SCC\nSBS4",
				   "Lymph-BNHL\nSBS9") 
								 
full_plot_dt$labels = paste0(full_plot_dt$cohort_name, 
							"\n", 
							full_plot_dt$signature)     
full_plot_dt$labels_to_keep = ifelse(full_plot_dt$labels %in% labels_to_keep, 
									 full_plot_dt$labels, 
									 "")           							 							 
pos = position_jitter(width = 0.3, seed = 5)                        
pdf(pff("data/004H_etiology_boxplots.pdf"), width = 6, height = 4)                        
ggboxplot(full_plot_dt, 
		  x = "etiologies", 
		  y = "Adj_R2", 
		  color = "etiologies", 
		  outlier.shape = NA)+
    geom_jitter(shape = 21, 
				aes(fill = etiologies, 
					size = mut_nums), 
				colour = "black", 
				position = pos)+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
         axis.title = element_text(size = 7),
         legend.text = element_text(size = 6),
         legend.title = element_text(size = 7))+              
    labs(x = "Signature etiology", 
		 y = "Adj. R2", 
		 size = "total mutations")+
    scale_colour_d3(drop = F)+
    scale_fill_d3(drop = F)+
    geom_text_repel(aes(label = labels_to_keep), 
					size = 2.5, 
					box.padding = 0.5, 
					min.segment.length = 0, 
					position = pos, 
					segment.size = 0.2)+
    scale_radius(range = c(1, 3))+
    guides(color = FALSE, fill = FALSE)+
	ylim(c(-0.1,1))
dev.off()                        

