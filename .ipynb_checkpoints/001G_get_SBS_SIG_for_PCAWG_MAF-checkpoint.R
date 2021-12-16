date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

library(stringi)

#Load in annotation file
annotation_file = fread("/.mounts/labs/reimandlab/private/generated_raw_data/PCAWG/PCAWG May2016 Data Release v1 %2F v1.1 %2F v1.2 %2F v1.3 %2F v1.4 - release_may2016v1.4 (1).csv")

#Load in PCAWG SBS probability file
PCAWG_sigs = fread("/.mounts/labs/reimandlab/private/generated_raw_data/PCAWG/SigProfilier_PCAWG_WGS_probabilities_SBS.csv")
PCAWG_sigs$Donor_ID = annotation_file$icgc_donor_id[match(PCAWG_sigs$Sample, 
														  annotation_file$tumor_wgs_icgc_specimen_id)]

#Load in PCAWG MAF
PCAWG_mutations_dt_hg38_SNP = fread(pff("/data/001C_PCAWG_mutations_hg38.csv"))[Variant_Type=="SNP"]

#Get trinucleotide context
PCAWG_mutations_dt_hg38_SNP$tri_nuc_context = toupper(substr(PCAWG_mutations_dt_hg38_SNP$ref_context, 10, 12))

#Get mutation type
substrRight = function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
PCAWG_mutations_dt_hg38_SNP$Mutation_type = substrRight(PCAWG_mutations_dt_hg38_SNP$Genome_Change, 3)

#Convert trinucleotide to pyramidine complement
where_reverse = substr(PCAWG_mutations_dt_hg38_SNP$tri_nuc_context, 2, 2) %in% c("A", "G")
PCAWG_mutations_dt_hg38_SNP$tri_nuc_context = ifelse(where_reverse, 
													 as.character(complement(DNAStringSet(PCAWG_mutations_dt_hg38_SNP$tri_nuc_context))), 
													 PCAWG_mutations_dt_hg38_SNP$tri_nuc_context)

#Convert mutation type to pyramidine complement
where_reverse = substr(PCAWG_mutations_dt_hg38_SNP$Mutation_type, 1, 1) %in% c("A", "G")
mutation_type_string = gsub(">","",PCAWG_mutations_dt_hg38_SNP$Mutation_type)
x = ifelse(where_reverse, 
		   as.character(complement(DNAStringSet(mutation_type_string))), 
		   mutation_type_string)
stri_sub(x, 2, 1) = ">" 
PCAWG_mutations_dt_hg38_SNP$Mutation_type = x

#Assign a mutational signature to every SNP (sig with most probability)
PCAWG_sigs = na.omit(PCAWG_sigs)
setkey(PCAWG_sigs,Donor_ID, `Mutation Type`, `Mutation Subtype`)

MATCHING_SIGS_DT = PCAWG_sigs[.(PCAWG_mutations_dt_hg38_SNP$Donor_ID, 
								PCAWG_mutations_dt_hg38_SNP$Mutation_type, 
								PCAWG_mutations_dt_hg38_SNP$tri_nuc_context),
							  .SD,.SDcols=-c(1,2,3,4,70)]

DT = data.table(value=unlist(MATCHING_SIGS_DT, 
							 use.names=FALSE), 
                colid = 1:nrow(MATCHING_SIGS_DT), 
				rowid = rep(names(MATCHING_SIGS_DT), 
							each=nrow(MATCHING_SIGS_DT)))
setkey(DT, colid, value)
t1 = DT[J(unique(colid)), rowid, mult="last"]

PCAWG_mutations_dt_hg38_SNP$SBS_Signature = t1

#Save signature file
fwrite(PCAWG_mutations_dt_hg38_SNP, pff("/data/001G_PCAWG_MAF_SNP_withSBSSIG.csv"))