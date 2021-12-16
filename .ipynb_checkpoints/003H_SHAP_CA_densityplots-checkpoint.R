source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

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
                                                 
plot_SHAP_plot = function(cohort_index, top_N_predictors){
   	
	cohort_name = colnames(importance)[-1][cohort_index]
	predictors = importance[[1]]
	
	#Remove hypermutated windows from lymphoma and leukemia
    if(cohort_name %in% c("Lymph-CLL", "Lymph-BNHL")){
        Preds = Preds[-c(292, 2438)]
    }
	
	#Get permutation test p-values for each predictor
    data = importance[[cohort_index + 1]]
    permutation_importances_dt = as.data.table(do.call("rbind.data.frame", 
									 lapply(list.files(p_val_paths[which(p_val_cancer_types == cohort_name)], 
													   full.name = T), fread)))
    predictor_pvals = unlist(lapply(1:length(predictors), 
									function(x)  sum(data[x] < permutation_importances_dt[[x]])/1000))
									
    significance = ifelse(predictor_pvals == 0, "*", "")
    names(significance) = predictors
	top_N_predictors = min(top_N_predictors, sum(significance == "*"))
	
	#Get predictor importance means and sd from bootstrap experiment for significant predictor                     
    Bootstrap_dt = as.data.table(do.call("rbind.data.frame", 
										 lapply(list.files(bootstrap_paths[which(bootstrap_cancer_types == cohort_name)], 
														   full.name = T), fread)))
    sig_predictor_means = unlist(lapply(which(significance == "*"), 
										function(x) mean(Bootstrap_dt[[x]])))                                   
    top_predictors = names(sig_predictor_means[order(sig_predictor_means, 
													 decreasing = T)][1:top_N_predictors])

	top_predictor_descriptions = predictor_descriptions$Predictor_descriptions[match(top_predictors, predictor_descriptions$Predictor_names)]	
    
    #Load in SHAP values
    SHAP_dt = fread(SHAP_paths[which(SHAP_cancer_types == cohort_name)])
										
    #Keep SHAP values for top predictors
    SHAP_dt_top = SHAP_dt[,.SD,.SDcols = top_predictors]
    colnames(SHAP_dt_top) = top_predictor_descriptions
    Input_top = Preds[,.SD,.SDcols = top_predictors]                                  
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
                                      

plot_dt = as.data.table(do.call("rbind.data.frame",lapply(seq(26)[-8], plot_SHAP_plot, 5)))

cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")					  

plot_dt_select = plot_dt[cohort_name %in% cancer_types_to_keep]

plot_dt_primarycancer = plot_dt_select[type == "Primary cancer"]
plot_dt_normal_tissue = plot_dt_select[type == "Normal tissue"]									
plot_dt_cancer_cell_line = plot_dt_select[type == "Cancer cell line"]
									
plot_dt_RT = plot_dt_select[type == "RT"]
phases = unlist(lapply(as.character(plot_dt_RT$variable), 
					   function(x) unlist(strsplit(x, split = " "))[2]))
plot_dt_RT_early = plot_dt_RT[which(phases %in% c("G1", "S1", "S2"))]
plot_dt_RT_late = plot_dt_RT[which(phases %in% c("S3", "S4", "G2"))]

									
#Create faceted plot
plot_dt_facet = as.data.table(rbind.data.frame(plot_dt_primarycancer,
											   plot_dt_normal_tissue,
											   plot_dt_cancer_cell_line,
											   plot_dt_RT_early,
											   plot_dt_RT_late))
									
plot_dt_facet$facet = c(rep("TCGA", nrow(plot_dt_TCGA)),
						rep("DNase-Seq", nrow(plot_dt_DNase)),
						rep("RT_Early", nrow(plot_dt_RT_early)), 
						rep("RT_Late", nrow(plot_dt_RT_late)))

									
pdf(pff("data/003H_SHAP_facet_densityplot.pdf"), height = 3, width = 8.5)                                      
ggplot(plot_dt_facet, aes(y = Input, x = SHAP))+
	facet_wrap(~type, nrow = 1)+
    geom_hex(bins = 50)+
    geom_smooth(method = "loess", se = F)+
    theme_bw()+
    scale_fill_gsea()+
	labs(y = "Input(log1p, scaled)")
dev.off()
									
									
									
pdf(pff("data/003H_CART_SHAP_densityplots/1.pdf"), height = 3, width = 3)                                      
ggplot(plot_dt_TCGA, aes(y = Input, x = SHAP))+
    geom_hex(bins = 50)+
    geom_smooth(method = "loess", se = F)+
    theme_bw()+
    scale_fill_gsea()+
	labs(y = "TCGA CA(log1p, scaled)")
dev.off()									
			
									
#Numbers
cor(plot_dt_TCGA$Input, plot_dt_TCGA$SHAP, method="spearman")
# [1] -0.7501813
cor.test(plot_dt_TCGA$Input, plot_dt_TCGA$SHAP, method="spearman")$p.val
# [1] 0
cor(plot_dt_DNase$Input, plot_dt_DNase$SHAP, method="spearman")
# [1] -0.7487669
cor.test(plot_dt_DNase$Input, plot_dt_DNase$SHAP, method="spearman")$p.val
# [1] 0	
cor(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method="spearman")
# [1] -0.7853214
cor.test(plot_dt_RT_early$Input, plot_dt_RT_early$SHAP, method="spearman")$p.val
# [1] 0									
cor(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method="spearman")
# [1] 0.7736708
cor.test(plot_dt_RT_late$Input, plot_dt_RT_late$SHAP, method="spearman")$p.val
# [1] 0													
									
									