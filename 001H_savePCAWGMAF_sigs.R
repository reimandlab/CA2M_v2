date_tag = "210317"
source(paste0("/.mounts/labs/reimandlab/private/users/oocsenas/CA2M/", 
# source(paste0("/.mounts/labs/reimandlab/private/users/jreimand/CA2M/", 
		date_tag, 
		"/bin/000_HEADER.R"))

source(pff("bin/999_process_data.R"))

library(stringi)

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

fwrite(PCAWG_mutations_dt_hg38_SNP, pff("data/001H_PCAWG_mutations_hg38_MAF_sigs.csv"))