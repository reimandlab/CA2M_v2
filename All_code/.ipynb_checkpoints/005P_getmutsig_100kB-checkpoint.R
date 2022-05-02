#Get cohort number as argument to loop through
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("000_HEADER.R")

source(pff("bin/999_process_data.R"))

#Load in PCAWG annotation file for sample codes
annotation_file = fread("")

#Load in COSMIC SBS probability file for all PCAWG mutations
PCAWG_sigs = fread("")

#Standardize donor ID's
PCAWG_sigs$Donor_ID = annotation_file$icgc_donor_id[match(PCAWG_sigs$Sample, 
														  annotation_file$tumor_wgs_icgc_specimen_id)]

#Load in PCAWG MAF file
PCAWG_mutations_dt_hg38_SNP = fread(pff("/data/001H_PCAWG_mutations_hg38_MAF_sigs.csv"))

#Get cohort name for cohort index argument
cancer_types_to_keep = c("Breast-AdenoCa", "Prost-AdenoCA", "Kidney-RCC", "Skin-Melanoma", 
						 "Uterus-AdenoCA","Eso-AdenoCa", 
						 "Stomach-AdenoCA","CNS-GBM", "Lung-SCC", "ColoRect-AdenoCA", "Biliary-AdenoCA", 
						 "Head-SCC", "Lymph-CLL", "Lung-AdenoCA",
						   "Lymph-BNHL",  "Liver-HCC", "Thy-AdenoCA")	
cohort_name = cancer_types_to_keep[cohort_index]

#Keep only mutations from that cohort
PCAWG_mutations_dt_hg38_SNP = PCAWG_mutations_dt_hg38_SNP[Project_Code == cohort_name]

#Get donors from that cohort
donors = unique(PCAWG_mutations_dt_hg38_SNP$Donor_ID)

#Create hg38 genome and bin into megabase-scale windows
genome = BSgenome.Hsapiens.UCSC.hg38
gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr", 1:22, sep = "")], 
							tilewidth = 100000, 
							cut.last.tile.in.chrom = TRUE)

#Get mutation counts for each genomic window for each signature in cohort of interest
#We use mut probability * number of mutations for each SBS to estimate this value
get_counts_signature = function(signature){ #Get counts for signature of interest

		print(signature)
	
		get_counts_sample = function(donor_ID_arg){ # Get counts for each donor in cohort
			
			get_counts_trinuc = function(tri_nuc_context_arg){ #Get counts for each trinucleotide context of donor

				get_counts_mut_type = function(mutaton_type_arg){ #Get counts for each mutation type of donor

					#Get probability of mutation for signature, donor, trinucleotide context, and mutation type
					prob = PCAWG_sigs[Donor_ID == donor_ID_arg & `Mutation Subtype` == tri_nuc_context_arg & `Mutation Type` == mutaton_type_arg][[signature]]
					
					#Get mutations for signature, donor, trinucleotide context, and mutation type
					muts = PCAWG_mutations_dt_hg38_SNP[Donor_ID == donor_ID_arg & tri_nuc_context == tri_nuc_context_arg & Mutation_type == mutaton_type_arg]
					
					#Get binned mutation counts signature, donor, trinucleotide context, and mutation type
					RMV = create_RMV(muts, 100000)[[3]]
					
					#Multiply binned mutation counts by probability of mutation to get estimated counts
					RMV_prob = RMV * prob

					return(RMV_prob)

				}
				
				#Combine genomic mutation counts for different mutation types for trinucleotide context
				counts_dt = do.call("cbind", lapply(unique(PCAWG_mutations_dt_hg38_SNP[Donor_ID == donor_ID_arg & tri_nuc_context == tri_nuc_context_arg]$Mutation_type), 
					get_counts_mut_type))
				
				#Sum estimated mutation counts for trinucleotide context
				counts = rowSums(counts_dt)

				return(counts)}
			
			#Combine genomic mutation counts for different trinucleotide contexts for donor
			counts_dt  = do.call("cbind", lapply(unique(PCAWG_mutations_dt_hg38_SNP[Donor_ID == donor_ID_arg]$tri_nuc_context), get_counts_trinuc))
			
			#Sum estimated mutation counts for trinucleotide context
			counts = rowSums(counts_dt)
			
			return(counts)}
		
		#Combine genomic mutation counts for donors for signature
		counts_dt  = do.call("cbind", lapply(donors, get_counts_sample))
		
		#Sum estimated mutation counts for signature
		counts = rowSums(counts_dt)
	
		#Round estimated mutation counts for signature
		counts_rounded = round(counts, digits = 0)
	
		return(counts_rounded)}

#Combine genomic mutation counts into one data table for different signatuers for cancer type
counts_dt = do.call("cbind", mclapply(colnames(PCAWG_sigs)[-c(1:4, 70)], get_counts_signature, mc.cores = 4))

#Add genomic coordinates
counts_dt_full = as.data.table(cbind.data.frame(as.character(seqnames(gr.windows)), 
																								start(gr.windows),
																								counts_dt))
colnames(counts_dt_full) = c("chr", "start", colnames(PCAWG_sigs)[-c(1:4, 70)])

#Keep only windows with >80% mappability
counts_dt_mappable = filter_windows_with_UMAP(counts_dt_full, 100000)

#Save genomic mutations counts for all signatures in cohort
fwrite(counts_dt_mappable, pff(paste0("data/005P_mutsigs_100KB/", cohort_name, ".csv")))

