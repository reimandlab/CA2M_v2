source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

#Load in required python packages to load and plot SHAP values using reticulate
library(reticulate)
shap = import("shap")
plt = import('matplotlib.pyplot')
backend=import("matplotlib.backends.backend_pdf")

#Load in top preds
importance = fread(pff("data/003A_feature_importances_allpreds.csv"))

#Load in predictor supplementary table
predictor_descriptions = fread(pff("data/001K_predictor_descriptions.csv"))           

#Load predictor matching list
matching_list = readRDS(pff("data/001J_Matching_predictors_dictionary.RDS"))

#Load in permutation test paths
p_val_paths = list.files(pff("data/003B_Feature_importance_permutations/") ,full.name = T)
p_val_cancer_types = unlist(lapply(list.files(pff("/data/003B_Feature_importance_permutations/")), 
											  function(x) unlist(strsplit(x, split = ".csv"))[1]))
                                 
#Load in bootstrap test paths                                 
bootstrap_paths = list.files(pff("data/003C_Feature_importance_bootstrapping/"), full.name = T)                
bootstrap_cancer_types = unlist(lapply(list.files(pff("data/003C_Feature_importance_bootstrapping/")), 
												  function(x) unlist(strsplit(x, split = ".csv"))[1]))
      
#Load in SHAP dataset paths                                     
SHAP_paths = list.files(pff("data/003F_RF_SHAP_Results/"), full.name = T)
SHAP_cancer_types = unlist(lapply(list.files(pff("data/003F_RF_SHAP_Results/")), 
								function(x) unlist(strsplit(x, split = ".csv"))[1]))

#Load in predictors
Preds = fread(pff("data/001G_All_preds_1MB.csv"))[, -c(1,2)]

								  
#Get data for SHAP vs. CA plot								  
plot_SHAP_plot = function(cohort_name){

	print(cohort_name)
	#Get predictor importances from permutation test
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
	#Get predictor names
	predictors = colnames(permutation_importances_dt)
	
	#Get predictor importances
    data = importance[[cohort_name]][match(predictors, importance[[1]])]
	
	#Get permutation test p-values for each predictor (fraction of times permuted importance is more significant than actual importance)
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))				
	significance = ifelse(predictor_pvals == 0, "Significant", "Non-Significant")
																

	#Remove hypermutated windows for Lymph CLL and BNHL
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
    }
			
    #Load in SHAP values
    SHAP_dt = fread(SHAP_paths[which(SHAP_cancer_types == cohort_name)])
			   
    Input_dt = as.data.table(apply(Preds, 2, function(x) scale(log(1+x)) ))
										
    #Create tall format                                  
    SHAP_dt.m = melt(SHAP_dt)
    Input_dt.m = melt(Input_dt)
                                      
    #Combine to final dt
    plot_dt.m = as.data.table(cbind.data.frame(variable = SHAP_dt.m$variable, 
											   Input = Input_dt.m$value, 
											   SHAP = SHAP_dt.m$value, 
                                             cohort_name = rep(cohort_name, 
																nrow(Input_dt.m)), 
											   significance = significance[match(SHAP_dt.m$variable, predictors)]))
								   
	
	#Downsample non-significant predictors						   
	set.seed(1)
	non_sig_to_keep = sample(predictors[which(predictor_pvals > 0.2)], sum(predictor_pvals == 0))															   
	plot_dt.m = plot_dt.m[variable %in% c(predictors[which(predictor_pvals == 0)], non_sig_to_keep)]							   
								   
								   
	plot_dt.m$variable = predictor_descriptions$Predictor_descriptions[match(plot_dt.m$variable, predictor_descriptions$Predictor_names)]
	plot_dt.m$type = predictor_descriptions$Predictor_categories[match(plot_dt.m$variable, predictor_descriptions$Predictor_descriptions)]
	
															   
								   
    return(plot_dt.m)}
                                      
#Get SHAP vs. CA data for all cancer types
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", 
						 "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", 
						 "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
plot_dt = as.data.table(do.call("rbind.data.frame",lapply(cancer_types_to_keep, plot_SHAP_plot)))

#Separate data based on predictor type 									
plot_dt_primarycancer = plot_dt[type == "Primary cancer"]
plot_dt_normal_tissue = plot_dt[type == "Normal tissue"]									
									
plot_dt_RT = plot_dt[type == "RT"]
phases = unlist(lapply(as.character(plot_dt_RT$variable), 
					   function(x) unlist(strsplit(x, split = " "))[2]))
plot_dt_RT_early = plot_dt_RT[which(phases %in% c("G1", "S1", "S2"))]
plot_dt_RT_late = plot_dt_RT[which(phases %in% c("S3", "S4", "G2"))]

									
#Create faceted plot for different predictor types (data.table long-form)
plot_dt_facet = as.data.table(rbind.data.frame(plot_dt_primarycancer,
											   plot_dt_normal_tissue,
											   plot_dt_RT_early,
											   plot_dt_RT_late))
									
plot_dt_facet$category = c(rep("Primary cancer", nrow(plot_dt_primarycancer)),
						rep("Normal tissue", nrow(plot_dt_normal_tissue)),
						rep("RT Early", nrow(plot_dt_RT_early)), 
						rep("RT Late", nrow(plot_dt_RT_late)))

					   
cor_dt = cbind.data.frame(category = c("Primary cancer", "Normal tissue", "RT Early", "RT Late", 
									   "Primary cancer", "Normal tissue", "RT Early", "RT Late"),
						  significance = c(rep("Significant", 4,), rep("Non-Significant", 4)),
						  label = c(cor(plot_dt_primarycancer[significance == "Significant"]$Input, 
										plot_dt_primarycancer[significance == "Significant"]$SHAP, method = "spearman"),
							cor(plot_dt_normal_tissue[significance == "Significant"]$Input, 
								plot_dt_normal_tissue[significance == "Significant"]$SHAP, method = "spearman"),
							cor(plot_dt_RT_early[significance == "Significant"]$Input, 
								plot_dt_RT_early[significance == "Significant"]$SHAP, method = "spearman"), 
							cor(plot_dt_RT_late[significance == "Significant"]$Input, 
								plot_dt_RT_late[significance == "Significant"]$SHAP, method = "spearman"), 
							cor(plot_dt_primarycancer[significance == "Non-Significant"]$Input, 
										plot_dt_primarycancer[significance == "Non-Significant"]$SHAP, method = "spearman"),
							cor(plot_dt_normal_tissue[significance == "Non-Significant"]$Input, 
								plot_dt_normal_tissue[significance == "Non-Significant"]$SHAP, method = "spearman"),
							cor(plot_dt_RT_early[significance == "Non-Significant"]$Input, 
								plot_dt_RT_early[significance == "Non-Significant"]$SHAP, method = "spearman"), 
							cor(plot_dt_RT_late[significance == "Non-Significant"]$Input, 
								plot_dt_RT_late[significance == "Non-Significant"]$SHAP, method = "spearman")))
					   
#Plot SHAP vs. CA density plot					   
plot_dt_facet$category = factor(plot_dt_facet$category, 
							 levels = c("Primary cancer", "Normal tissue", "RT Early", "RT Late"))		
plot_dt_facet$significance = factor(plot_dt_facet$significance, 
							 levels = c("Significant", "Non-Significant"))	
					   
cor_dt$category = factor(cor_dt$category, 
							 levels = c("Primary cancer", "Normal tissue", "RT Early", "RT Late"))		
cor_dt$significance = factor(cor_dt$significance, 
							 levels = c("Significant", "Non-Significant"))	
					   
pdf(pff("data/003I_SHAP_facet_densityplot_nonsignificant.pdf"))                                      
ggplot((plot_dt_facet), aes(y = Input, x = SHAP))+
	facet_grid(significance ~ category, scales = "free", switch = "y")+ #Facet based on predictor type
    geom_hex(bins = 50)+ #Add hexagonal density plot
    geom_smooth(method = "loess",span = 0.6, se = FALSE)+ #Add loess smooth line
    theme_bw()+
    scale_fill_gsea()+ #Control colours
	labs(y = "Input(log1p, scaled)") + 
	geom_text(data = cor_dt,
		  mapping = aes(x = -Inf, y = Inf, label = paste0("rho=", round(label,2))),
		  hjust   = -0.1,
		  vjust   = 1.4,
			color = "red")
dev.off()

#Correlations for the different plots
cor(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "spearman")
# [1] -0.1079348
cor.test(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "spearman")$p.val
# [1] 0
cor(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "spearman")
# [1] -0.03155112
cor.test(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "spearman")$p.val
# [1] 0	
cor(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "spearman")
# [1] -0.7053419
cor.test(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "spearman")$p.val
# [1] 0									
cor(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "spearman")
# [1] 0.4395361
cor.test(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "spearman")$p.val
# [1] 0													
									
									