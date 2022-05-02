source("000_HEADER.R")

source(pff("bin/999_process_data.R"))

#Load in PCAWG MAF table, hg38 and keep SNV's only
PCAWG_mutations_dt_hg38_SNP = fread(pff("/data/001C_PCAWG_mutations_hg38.csv"))[Variant_Type=="SNP"]

#Get trinucleotide context for each SNV
PCAWG_mutations_dt_hg38_SNP$tri_nuc_context = toupper(substr(PCAWG_mutations_dt_hg38_SNP$ref_context, 10, 12))

#Get mutation type for each SNV
substrRight = function(x, n){ #Function takes last n characters of any string
  substr(x, nchar(x)-n+1, nchar(x))
}
PCAWG_mutations_dt_hg38_SNP$Mutation_type = substrRight(PCAWG_mutations_dt_hg38_SNP$Genome_Change, 3)

#Convert all trinucleotides to pyrimidine complement
where_reverse = substr(PCAWG_mutations_dt_hg38_SNP$tri_nuc_context, 2, 2) %in% c("A", "G") #Get purine complements
PCAWG_mutations_dt_hg38_SNP$tri_nuc_context = ifelse(where_reverse,  #Convert these to pyrimidine complement
													 as.character(complement(DNAStringSet(PCAWG_mutations_dt_hg38_SNP$tri_nuc_context))), 
													 PCAWG_mutations_dt_hg38_SNP$tri_nuc_context)

#Convert all mutation types to pyramidine complement
where_reverse = substr(PCAWG_mutations_dt_hg38_SNP$Mutation_type, 1, 1) %in% c("A", "G") #Get purine complements
mutation_type_string = gsub(">", "", PCAWG_mutations_dt_hg38_SNP$Mutation_type)
x = ifelse(where_reverse, #Convert these to pyrimidine complement
		   as.character(complement(DNAStringSet(mutation_type_string))), 
		   mutation_type_string)
stri_sub(x, 2, 1) = ">" 
PCAWG_mutations_dt_hg38_SNP$Mutation_type = x

#Save PCAWG MAF file with trinucleotide context and mutation type columns
fwrite(PCAWG_mutations_dt_hg38_SNP, pff("data/001H_PCAWG_mutations_hg38_MAF_sigs.csv"))