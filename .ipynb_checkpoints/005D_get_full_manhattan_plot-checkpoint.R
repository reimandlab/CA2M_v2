date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

input_data_dir = "/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/INPUT_DATA/"

#Load in error dt paths
#Get dt for 12 core cancer types
cancer_types_to_keep = c("Breast-AdenoCa", "Kidney-RCC", "CNS-GBM", "Eso-AdenoCa",
						 "Stomach-AdenoCA","Lung-SCC", 
						 "ColoRect-AdenoCA",  
						  "Biliary-AdenoCA", "Head-SCC" ,"Lung-AdenoCA", 
						 "Thy-AdenoCA", "Liver-HCC", "Skin-Melanoma")	
q_val_dt = fread(pff("/data/005C_qvaldt_100KBerrorwindows.csv"))[,.SD,.SDcols = c("chr", "start", cancer_types_to_keep)]

#Convert chr to numeric
q_val_dt$chr = as.numeric(unlist(lapply(q_val_dt$chr, 
										function(x) (unlist(strsplit(x, split = "chr"))[2]))))

#Melt data table
q_val_dt.m = melt(q_val_dt, id.vars = c("chr", "start"))

#Cap significance for visualization purposes
q_val_dt.m$value[which(q_val_dt.m$value < 10^-20)] = 10^-20                                      
                                      
#Get significant points
q_val_dt.m$sig = ifelse(q_val_dt.m$value < 0.05, "yes", "no")                                      

#Get PCAWG colours
PCAWG_colours = readRDS(paste0(input_data_dir, "PCAWG_colour_palette.RDS"))
q_val_dt.m$colour_sig = ifelse(q_val_dt.m$sig == "yes", 
							   as.character(q_val_dt.m$variable), 
							   "gray")                                      
                                      
#Create cumulative base pair position
q_val_dt.m$cum_bp = q_val_dt.m$start + hg38$cumlen[match(q_val_dt.m$chr, 
														 hg38$chrom)]       

# q_val_dt.m = na.omit(q_val_dt.m)
q_val_dt.m$colour_sig = factor(q_val_dt.m$colour_sig)                                
q_val_dt.m$sig = factor(q_val_dt.m$sig)
                                      
#PCAWG palette                                      
PCAWG_pal = as.character(PCAWG_colours)[match(tolower(levels(q_val_dt.m$colour_sig)), 
											  names(PCAWG_colours))]
PCAWG_pal[6] =  "#D3D3D3"
                                  
p1_rast = ggplot(q_val_dt.m, aes(x = cum_bp, y = -1*log10(value)))+
    geom_point_rast(data = q_val_dt.m[sig=="no"], 
					shape = 21, 
					size = 1.5, 
					alpha = 0.75, 
					aes(fill = colour_sig, colour = sig))+
    geom_point(data = q_val_dt.m[sig=="yes"], 
			   shape = 21,
			   size = 1.5,
			   alpha = 0.75,
			   aes(fill = colour_sig, colour = sig))+
    theme_bw()+
    scale_fill_manual(values = PCAWG_pal, 
					  breaks = levels(q_val_dt.m$colour_sig)[-6],
					  guide = F)+
    scale_colour_manual(values = c("grey", "black"), 
						guide = F)+
    geom_hline(yintercept = -1*log10(0.05), 
			   colour = "red",
			   linetype = "dashed")+
    scale_x_continuous(breaks = (hg38$cumlen[1:22] + hg38$cumlen[2:23])/2, 
					   labels = 1:22)+
    labs(x = "Chromosome",
		 fill = "Project Code",
		 y = "-log10 (q-value)")+
    theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(vjust = 8, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9),    
    axis.ticks.x = element_blank())+
    geom_vline(xintercept = hg38$cumlen[1:23],
			   linetype = "dashed",
			   colour = "grey",
			   alpha = 0.5)
                                      
pdf(pff("/data/005D_14cancertypes_manhattan_plot_sig_rast.pdf"), width = 9.5, height = 4)                             
p1_rast                                      
dev.off()   
                                                                          
#Get legend
p1_legend = ggplot(q_val_dt.m, aes(x = cum_bp, y = -1*log10(value)))+
    geom_point_rast(data = q_val_dt.m[sig=="no"], 
					shape = 21, 
					size = 1.5, 
					alpha = 0.75, 
					aes(fill = colour_sig, colour = sig))+
    geom_point(data = q_val_dt.m[sig=="yes"], 
			   shape = 21,
			   size = 1.5,
			   alpha = 0.75,
			   aes(fill = colour_sig, colour = sig))+
    theme_bw()+
    scale_fill_manual(values = PCAWG_pal, 
					  breaks = levels(q_val_dt.m$colour_sig)[-6])+
    scale_colour_manual(values = c("grey", "black"), 
						guide = F)+
    geom_hline(yintercept = -1*log10(0.05), 
			   colour = "red",
			   linetype = "dashed")+
    scale_x_continuous(breaks = (hg38$cumlen[1:22] + hg38$cumlen[2:23])/2, 
					   labels = 1:22)+
    labs(x = "Chromosome",
		 fill = "Project Code",
		 y = "-log10 (q-value)")+
    theme( 
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(vjust = 8, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9),    
    axis.ticks.x = element_blank())+
    geom_vline(xintercept = hg38$cumlen[1:23],
			   linetype = "dashed",
			   colour = "grey",
			   alpha = 0.5)										
										
legend_full <- cowplot::get_legend(p1_legend)

pdf(pff("/data/005D_14cancertypes_manhattan_plot_LEGEND.pdf"), width = 4, height = 4)
grid.newpage()             
grid.draw(legend_full)                                      
dev.off()                                      
