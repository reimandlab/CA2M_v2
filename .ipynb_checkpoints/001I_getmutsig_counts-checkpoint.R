#Get cohort number as argument
args = commandArgs(trailingOnly = TRUE)
cohort_index = as.integer(args[1])

source("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M_v2/bin/000_HEADER.R")

source(pff("bin/999_process_data.R"))

library(stringi)

#Load in annotation file
annotation_file = fread("/.mounts/labs/reimandlab/private/generated_raw_data/PCAWG/PCAWG May2016 Data Release v1 %2F v1.1 %2F v1.2 %2F v1.3 %2F v1.4 - release_may2016v1.4 (1).csv")

#Load in PCAWG SBS probability file
PCAWG_sigs = fread("/.mounts/labs/reimandlab/private/generated_raw_data/PCAWG/SigProfilier_PCAWG_WGS_probabilities_SBS.csv")
PCAWG_sigs$Donor_ID = annotation_file$icgc_donor_id[match(PCAWG_sigs$Sample, 
														  annotation_file$tumor_wgs_icgc_specimen_id)]

#Load in PCAWG MAF
PCAWG_mutations_dt_hg38_SNP = fread(pff("/data/001H_PCAWG_mutations_hg38_MAF_sigs.csv"))
cohort_name = unique(PCAWG_mutations_dt_hg38_SNP$Project_Code)[cohort_index]
PCAWG_mutations_dt_hg38_SNP = PCAWG_mutations_dt_hg38_SNP[Project_Code == cohort_name]
donors = unique(PCAWG_mutations_dt_hg38_SNP$Donor_ID)

#Get counts for each cohort
genome = BSgenome.Hsapiens.UCSC.hg38
gr.windows = tileGenome(seqinfo(Hsapiens)[paste("chr", 1:22, sep = "")], 
							tilewidth = 1000000, 
							cut.last.tile.in.chrom = TRUE)


		
get_counts_signature = function(signature){

		print(signature)
		get_counts_sample = function(donor_ID_arg){
			

			get_counts_trinuc = function(tri_nuc_context_arg){

				get_counts_mut_type = function(mutaton_type_arg){

					#Get probability of mutation for signature, donor, and mutation type
					prob = PCAWG_sigs[Donor_ID == donor_ID_arg & `Mutation Subtype` == tri_nuc_context_arg & `Mutation Type` == mutaton_type_arg][[signature]]

					muts = PCAWG_mutations_dt_hg38_SNP[Donor_ID == donor_ID_arg & tri_nuc_context == tri_nuc_context_arg & Mutation_type == mutaton_type_arg]

					RMV = create_RMV(muts, 1000000)[[3]]

					RMV_prob = RMV * prob

					return(RMV_prob)

				}

				counts_dt = do.call("cbind", lapply(unique(PCAWG_mutations_dt_hg38_SNP[Donor_ID == donor_ID_arg & tri_nuc_context == tri_nuc_context_arg]$Mutation_type), 
					get_counts_mut_type))

				counts = rowSums(counts_dt)

				return(counts)}

			counts_dt  = do.call("cbind", lapply(unique(PCAWG_mutations_dt_hg38_SNP[Donor_ID == donor_ID_arg]$tri_nuc_context), get_counts_trinuc))
			counts = rowSums(counts_dt)
			return(counts)}
		
		counts_dt  = do.call("cbind", lapply(donors, get_counts_sample))
		counts = rowSums(counts_dt)

		counts_rounded = round(counts, digits = 0)
		return(counts_rounded)}
counts_dt = do.call("cbind", mclapply(colnames(PCAWG_sigs)[-c(1:4, 70)], get_counts_signature, mc.cores = 4))
counts_dt_full = as.data.table(cbind.data.frame(as.character(seqnames(gr.windows)), 
																								start(gr.windows),
																								counts_dt))
colnames(counts_dt_full) = c("chr", "start", colnames(PCAWG_sigs)[-c(1:4, 70)])


counts_dt_mappable = filter_windows_with_UMAP(counts_dt_full, 1000000)

fwrite(counts_dt_mappable, pff(paste0("data/001I_PCAWG_sigs_new/", cohort_name, ".csv")))

