source("000_HEADER.R")

source(pff("/bin/999_run_randomforest_experiment.R"))

#Get normal tissue CA + RT Predictors
Normal_preds = fread(pff("data/001G_Normal_preds_1MB.csv"))
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]
Normal_preds = as.data.table(cbind.data.frame(Normal_preds[,-c(1,2)], RT))

#Get primary cancer CA + RT Predictors
Tumor_preds = fread(pff("data/001G_PrimaryTumor_preds_1MB.csv"))
RT = fread(pff("data/001F_ENCODE_repliseq_1MBwindow_processed.csv"))[,-c(1,2)]
Tumor_preds = as.data.table(cbind.data.frame(Tumor_preds[,-c(1,2)], RT))

#Load in data paths for observed vs. expected RF outputs
tumour_paths = list.files(pff("data/002D_tumourCAplusRT_RF_obsvsexpected"), full.name = T)
normal_paths = list.files(pff("data/002E_normalCAplusRT_RF_obsvsexpected"), full.name = T)

#Get cancer types of interest for plotting
cohort_names = unlist(lapply(list.files(pff("data/002D_tumourCAplusRT_RF_obsvsexpected")),
                           function(x) unlist(strsplit(x, split = ".csv"))[1]))

cohorts_to_examine = c("PANCAN", "Breast-AdenoCa", "Prost-AdenoCA", "Eso-AdenoCa", "Liver-HCC")
                           
#Load in observed vs. expected values and put into data table for cancer types of interest
get_dt = function(cohort_name, paths, Preds, type){
        
    cohort_index = which(cohort_names == cohort_name) #Get index of cohort of interest
    
    data = fread(paths[cohort_index]) #Load in observed vs. expected
    
    R2 = (cor(data$observed,data$predicted))**2 #Get R2 accuracy score
    
    adjR2 = adjust_R2(R2, nrow(Preds), ncol(Preds)) #Get adjusted R2 accuracy score
    
	#Combine required data into one data table
    dt = as.data.table(cbind.data.frame(observed = data$observed, 
										predicted = data$predicted, 
										cancer_type = cohort_name, 
										predictor_set = paste0(cohort_name, 
															  "\n", type, " CA + RT\n", 
															  "Adj.R2=", 
															  round(adjR2, 2))))
    
    return(dt)

}

#Get required data for models trained on primary cancer CA + RT
tumour_dt = as.data.table(do.call("rbind.data.frame", 
								  lapply(cohorts_to_examine, get_dt, tumour_paths, Tumor_preds, "Tumour")))

#Get required data for models trained on normal tissue CA + RT
normal_dt = as.data.table(do.call("rbind.data.frame", 
								  lapply(cohorts_to_examine, get_dt, normal_paths, Normal_preds, "Normal")))
							 
#Combine into one data table                          
full_dt = as.data.table(rbind.data.frame(tumour_dt, normal_dt))                           
                           
                           
#Plot facet wrapped scatter plots
full_dt$predictor_set = factor(full_dt$predictor_set, 
							   levels = unique(full_dt$predictor_set))

###################################Functions to change scales on facets                           
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}
                           
CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)
                           
facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}                           

################################################################################################################                           
                           
#Create observed vs. expected scatterplots for cancer/normal tissue RF models                           
pdf(pff("data/002F_facetwrap_errorplots.pdf"), width = 5.5, height = 3)
ggplot(full_dt, aes(x = observed, y = predicted))+
    rasterise(geom_point(alpha = 0.4, size = 0.3), dpi = 400)+ #Rasterize scatterplots
    geom_smooth(method = 'loess', span = 0.9, size = 0.5)+ #Add loess line                     
    theme_bw()+
    labs(x = "Observed mutations per Mb", y = "Predicted mutations per Mb")+ #Add labels
    theme(axis.title = element_text(size = 7), #Control text
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, colour = "black", angle = 30, hjust = 1),
         strip.text = element_text(size = 6))+
    facet_wrap(~predictor_set, nrow = 2, scales = "free")+ #Facet based on predictors and cancer type                     
    facet_wrap_custom(~predictor_set, nrow = 2, scales = "free", scale_overrides = list(  #Control x and y scales for each facet                      
        scale_override(1, scale_y_continuous(breaks = seq(4000, 20000, 4000), limits=c(4000, 20000))),
        scale_override(2, scale_y_continuous(breaks = seq(300, 900, 200), limits=c(300, 950))),
        scale_override(3, scale_y_continuous(breaks = seq(100, 600, 100), limits=c(100, 600))),
        scale_override(4, scale_y_continuous(breaks = seq(0, 3000, 1000), limits=c(0, 3000))),                      
        scale_override(5, scale_y_continuous(breaks = seq(500, 3000, 500), limits=c(500, 3000))),                     
        scale_override(6, scale_y_continuous(breaks = seq(4000, 20000, 4000), limits=c(4000, 20000))),
        scale_override(7, scale_y_continuous(breaks = seq(300, 900, 200), limits=c(300, 950))),
        scale_override(8, scale_y_continuous(breaks = seq(100, 600, 100), limits=c(100, 600))),
        scale_override(9, scale_y_continuous(breaks = seq(0, 3000, 1000), limits=c(0, 3000))),                      
        scale_override(10, scale_y_continuous(breaks = seq(500, 3000, 500), limits=c(500, 3000)))))                          
dev.off() 