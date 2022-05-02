source("000_HEADER.R")

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
    
	#Rank significant predictors by bootstrap mean
	top_predictors = predictors[which(predictor_pvals == 0)]
					
	#Get predictor descriptions								   
	top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_predictors, predictor_descriptions$Predictor_names)]	
	
    #Load in SHAP values
    SHAP_dt = fread(SHAP_paths[which(SHAP_cancer_types == cohort_name)])
	
	#Remove hypermutated windows for Lymph CLL and BNHL
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
    }
									   
    #Keep SHAP values for top predictors
    SHAP_dt_top = SHAP_dt[,.SD,.SDcols = top_predictors]
    colnames(SHAP_dt_top) = top_predictor_descriptions
    Input_top = as.data.table(apply(Preds[,.SD, .SDcols = top_predictors], 2, function(x) scale(log(1+x)) ))
    colnames(Input_top) = top_predictor_descriptions  
										
    #Create tall format                                  
    SHAP_dt.m = melt(SHAP_dt_top)
    Input_dt.m = melt(Input_top)
                                      
    #Combine to final dt
    plot_dt.m = as.data.table(cbind.data.frame(variable = SHAP_dt.m$variable, 
											   Input = Input_dt.m$value, 
											   SHAP = SHAP_dt.m$value, 
                                             cohort_name = rep(cohort_name, 
																nrow(Input_dt.m))))
									
	plot_dt.m$type = predictor_descriptions$Predictor_categories[match(plot_dt.m$variable, predictor_descriptions$Predictor_descriptions)]									
    
    return(plot_dt.m)}
                                      
#Get SHAP vs. CA data for all cancer types
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
plot_dt = as.data.table(do.call("rbind.data.frame",lapply(cancer_types_to_keep, plot_SHAP_plot)))

#Separate data based on predictor type 									
plot_dt_primarycancer = plot_dt[type == "Primary cancer"]
plot_dt_normal_tissue = plot_dt[type == "Normal tissue"]									
plot_dt_cancer_cell_line = plot_dt[type == "Cancer cell line"]
									
plot_dt_RT = plot_dt[type == "RT"]
phases = unlist(lapply(as.character(plot_dt_RT$variable), 
					   function(x) unlist(strsplit(x, split = " "))[2]))
plot_dt_RT_early = plot_dt_RT[which(phases %in% c("G1", "S1", "S2"))]
plot_dt_RT_late = plot_dt_RT[which(phases %in% c("S3", "S4", "G2"))]

									
#Create faceted plot for different predictor types (data.table long-form)
plot_dt_facet = as.data.table(rbind.data.frame(plot_dt_primarycancer,
											   plot_dt_normal_tissue,
											   plot_dt_cancer_cell_line,
											   plot_dt_RT_early,
											   plot_dt_RT_late))
									
plot_dt_facet$facet = c(rep("Primary cancer", nrow(plot_dt_primarycancer)),
						rep("Normal tissue", nrow(plot_dt_normal_tissue)),
						rep("Cancer cell line", nrow(plot_dt_cancer_cell_line)),
						rep("RT Early", nrow(plot_dt_RT_early)), 
						rep("RT Late", nrow(plot_dt_RT_late)))
					   
cor_dt = cbind.data.frame(facet = c("Primary cancer", "Normal tissue", "RT Early", "RT Late"),
						  label = c(cor(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "pearson"),
							cor(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "pearson"),
							cor(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "pearson"), 
							cor(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "pearson")))

#Plot SHAP vs. CA density plot					   
plot_dt_facet$facet = factor(plot_dt_facet$facet, 
							 levels = c("Primary cancer", "Normal tissue", "RT Early", "RT Late"))									
pdf(pff("data/003H_SHAP_facet_densityplot.pdf"), height = 3, width = 8.5)                                      
ggplot(na.omit(plot_dt_facet), aes(y = Input, x = SHAP))+
	facet_wrap(~facet, scales = "free", nrow = 1, drop = T)+ #Facet based on predictor type
    geom_hex(bins = 50)+ #Add hexagonal desnity plot
    geom_smooth(method = "lm", se = F)+ #Add loess smooth line
    theme_bw()+
    scale_fill_gsea()+ #Control colours
	labs(y = "Input(log1p, scaled)") + 
	geom_text(data = cor_dt,
		  mapping = aes(x = -Inf, y = Inf, label = paste0("r=", round(label,2))),
		  hjust   = -0.1,
		  vjust   = 1.4,
			color = "red")
dev.off()

#Correlations for the different plots
cor(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "spearman")
# [1] -0.6626878
cor.test(plot_dt_primarycancer$Input, plot_dt_primarycancer$SHAP, method = "spearman")$p.val
# [1] 0
cor(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "spearman")
# [1] -0.7203573
cor.test(plot_dt_normal_tissue$Input, plot_dt_normal_tissue$SHAP, method = "spearman")$p.val
# [1] 0	
cor(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "spearman")
# [1] -0.7807336
cor.test(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method = "spearman")$p.val
# [1] 0									
cor(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "spearman")
# [1] 0.7109059
cor.test(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method = "spearman")$p.val
# [1] 0													
									
									