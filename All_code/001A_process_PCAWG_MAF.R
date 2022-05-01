source("000_HEADER.R")

#Load in PCAWG MAF file
PCAWG_mutations_dt = fread(PCAWG_MAF_path, fill = T, showProgress = T)

#Keep only SNV's
PCAWG_mutations_dt_SNV = PCAWG_mutations_dt[Variant_Type == "SNP"]

#Remove hypermutated samples (>90k mutations)
mutations_per_sample = table(PCAWG_mutations_dt_SNV$Donor_ID)
hypermutated_samples = names(which(mutations_per_sample > 90000))
PCAWG_mutations_dt_SNV_nohyper = PCAWG_mutations_dt_SNV[-which(PCAWG_mutations_dt_SNV$Donor_ID %in% hypermutated_samples)]

#Normalize chromosome number/name
seqnames_PCAWG = paste("chr", PCAWG_mutations_dt_SNV_nohyper$Chromosome, sep = "")

#Create GRanges object for hg19 mutations
mutation_ranges_hg19 = GRanges(seqnames = seqnames_PCAWG, 
							   IRanges(PCAWG_mutations_dt_SNV_nohyper$Start_position, 
									   PCAWG_mutations_dt_SNV_nohyper$End_position), 
							   index = 1:nrow(PCAWG_mutations_dt_SNV_nohyper))

#Import liftover files
path = system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(paste0(input_data_dir, "/hg19ToHg38.over.chain"))

##Liftover mutations to hg38 from hg19
mutation_ranges_hg38 = unlist(liftOver(mutation_ranges_hg19, ch))

#Convert genomic coordinates in MAF file to hg38
PCAWG_mutations_dt_hg38 = PCAWG_mutations_dt_SNV_nohyper[mutation_ranges_hg38$index]
PCAWG_mutations_dt_hg38$Chromosome = as.character(seqnames(mutation_ranges_hg38))
PCAWG_mutations_dt_hg38$Start_position = start(mutation_ranges_hg38)
PCAWG_mutations_dt_hg38$End_position = end(mutation_ranges_hg38)

#Save hg38 PCAWG MAF file
fwrite(PCAWG_mutations_dt_hg38, pff("/data/001A_PCAWG_mutations_hg38.csv"))