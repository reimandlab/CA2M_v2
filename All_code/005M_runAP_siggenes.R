source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/INPUT_DATA/"

gmt = paste0(input_data_dir, "GOBP_REAC.gmt")

#Load in gene pval dt and remove genes with NAs
gene_pval_dt = na.omit(fread(pff("data/005L_gene_pvaldt_100KBerrorwindows.csv")))


cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	

gene_pval_dt_core = gene_pval_dt[, .SD, .SDcols = cancer_types_to_keep]

gene_pval_dt_core_matrix = data.matrix(gene_pval_dt_core)
rownames(gene_pval_dt_core_matrix) =  gene_pval_dt[[1]]

AP_integrated = ActivePathways(unique(gene_pval_dt_core_matrix), 
							 gmt, 
							 correction.method="fdr", 
							 cytoscape.file.tag = pff("data/005M_Integrated_AP_cytofiles/"))

fwrite(unique(gene_pval_dt_core_matrix), pff("data/005G_intergratedAP_genepval_table.csv"))
fwrite(AP_integrated, pff("data/005G_integratedAP_results.csv"))

nrow(AP_integrated)
# [1] 177

#Get how many cancer types each pathway takes evidence from
num_evidence = sapply(1:nrow(AP_integrated), 
					  function(x) length(unlist(strsplit(unlist(AP_integrated$evidence[x]), 
														 split = ",", fixed = T))))
num_combined = sum(unlist(AP_integrated$evidence) == "combined")
					  
sum(num_evidence > 1) + num_combined				  
# [1] 142
					  
#Edit subroups file to add PCAWG colours
library(mgsub)

subgroups = fread(pff("data/005M_Integrated_AP_cytofiles/subgroups.txt"))

change_cols = function(file, old_cols, new_cols){
    
    file$instruct = unlist(lapply(file$instruct, function(x) mgsub(x, old_cols, new_cols)))
                                
    return(file)                            
    
    
}
PCAWG_colours = readRDS(paste0(input_data_dir, "PCAWG_colour_palette.RDS"))

new_cols = c(as.character(PCAWG_colours)[match(tolower(colnames(subgroups)[-c(1,  
																			  ncol(subgroups))]), 
											   names(PCAWG_colours))])
new_cols[12] = "#FF0000"								  
								  
old_cols = c('#FF0000','#FF6000','#FFBF00','#DFFF00','#80FF00','#20FF00',
			 '#00FF40','#00FF9F','#00FFFF','#009FFF','#0040FF','#2000FF',
			 '#8000FF','#DF00FF','#FF00BF','#FF0060')
                                
subgroups_PCAWG = change_cols(subgroups, old_cols, new_cols)
write.table(subgroups_PCAWG, 
			pff("data/005M_Integrated_AP_cytofiles/subgroups_PCAWG.txt"), 
			quote = FALSE, 
			sep = '\t', 
		    row.names = F) 
								  
#Create legend for enrichment map
color_table = as.data.table(cbind.data.frame(group = colnames(subgroups)[-c(1,ncol(subgroups))], color = new_cols))
color_table = color_table[order(group)]								  
								  
dummy_plot = ggplot(data.frame(color_table,
							 "value" = 1), 
					aes(x = group, y = value, fill = group)) +
	geom_bar(stat = "identity") +
	scale_fill_manual(name = "Contribution", breaks = color_table$group, values = color_table$color)

legend_full <- cowplot::get_legend(dummy_plot)
pdf(pff("/data/005M_Integrated_AP_cytofiles/legend_PCAWG.pdf"))
grid.newpage()             
grid.draw(legend_full)                                      
dev.off()                                      

