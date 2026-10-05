# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_sensitivity_analysis_mom
# Description: Sensitivity analyses #1 and #2 on the female-specific bmi-associated SNPs
# ------------------------------------------------------------------

# --------------------------
# Load Libraries
# --------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rtracklayer")
library(rtracklayer)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(reshape2)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)


# Female BMI-associated SNPs (columns: 'variant','CHR', "POS",'Tested_Allele','BETA')
female_bmi_snps <- freadf('exposure_file.txt')

# --------------------------
# Filter 1
# --------------------------
freadf <- function(df){return(as.data.frame(fread(df)))} # Helper function: Load file as data.frame
fg_outcome <- freadf("FG_EUR_GCST90002232.tsv")
fi_outcome <- freadf("FI_EUR_GCST90002238.tsv")
hba1c_outcome <- freadf("HbA1c_EUR_GCST90002244.tsv")
homair_outcome <- freadf("HOMAIR_GCST90267567.tsv")
# homair_outcome is in GRCh38, the rest is in GRCh37

# Find significant snps associated with outcome (5e-8 filter)
filter_snps <- function(df) {
  df$p_value <- as.numeric(df$p_value) # Make sure p-value is numeric 
  df[df$p_value < 5e-8, ] # genome-wide filtering
}

fg_outcome <- filter_snps(fg_outcome)
fi_outcome <- filter_snps(fi_outcome)
hba1c_outcome <- filter_snps(hba1c_outcome)
homair_outcome <- filter_snps(homair_outcome)

# Create CHR:POS variant column
fg_outcome$variant <- paste(fg_outcome$chromosome, fg_outcome$base_pair_location, sep=":")
fi_outcome$variant <- paste(fi_outcome$chromosome, fi_outcome$base_pair_location, sep=":")
hba1c_outcome$variant <- paste(hba1c_outcome$chromosome, hba1c_outcome$base_pair_location, sep = ":")
homair_outcome$variant <- paste(homair_outcome$chromosome, homair_outcome$base_pair_location, sep = ":")

# Original PRS from exposure data 
PRS_original <- female_bmi_snps[,c('variant','CHR', "POS",'Tested_Allele','BETA')]

# Removing the snps from original PRS
PRS_outcome <- PRS_original[!PRS_original$variant %in% fg_outcome$variant, ]
PRS_outcome <- PRS_outcome[!PRS_outcome$variant %in% fi_outcome$variant, ]
PRS_outcome <- PRS_outcome[!PRS_outcome$variant %in% hba1c_outcome$variant, ]

# ---------------------
# Liftover of PRS variant coordinates from GRCh37 to GRCh38
# ---------------------

# Import chain file for GRCh37 to GRCh38 conversion
chain <- import.chain("hg19ToHg38.over.chain")

# Add "chr" prefix to chromosomes (required for GRanges)
PRS_outcome$CHR <- ifelse(grepl("^chr",PRS_outcome$CHR),
                        PRS_outcome$CHR,
                        paste0("chr",PRS_outcome$CHR))

# Ensure base pair positions are numeric
PRS_outcome$POS <- as.numeric(as.character(PRS_outcome$POS))

# Create a GRanges object for liftover
gr <- GRanges(
  seqnames = PRS_outcome$CHR,
  ranges = IRanges(start = PRS_outcome$POS, end = PRS_outcome$POS),
  chromosome = PRS_outcome$CHR, 
  base_pair_location = PRS_outcome$POS, 
  effect_allele = PRS_outcome$Tested_Allele, 
  beta = PRS_outcome$BETA, 
)
print(head(gr))

lifted <- liftOver(gr, chain) # Perform conversion
lifted_gr <- unlist(lifted) # Flatten result into a GRanges object
lifted <- as.data.frame(lifted_gr)

# Remove 'chr' prefix for compatibility with downstream datasets
lifted$chromosome <- sub("^chr", "", lifted$chromosome) # Removing the 'chr' prefix from chromosome number

# Create a new variant ID in the format CHR:POS
lifted$variant_id <- paste0(lifted$chromosome, ":", lifted$base_pair_location) 

# Remove overlapping SNPs with HOMA-IR 
lifted <- lifted[!lifted$variant_id %in% homair_outcome$variant, ] # No SNPs removed

# Save output file - Model 2
write.table(PRS_outcome[, c("variant", "Tested_Allele", "BETA")], file = "bmi_snp_mom_model2.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# --------------------------
# Filter 2 
# --------------------------
child_BMI <- read.csv("bmi_snp_child_model1.txt", sep = "\t") # 34 child BMI-associated SNPs
child_BMI$variant <- child_BMI$V1 # first column 

PRS_outcome <- PRS_outcome[!PRS_outcome$variant %in% child_BMI$variant, ] # 4 SNPs removed 

# Save output file - Model 3
write.table(PRS_outcome[, c("variant", "Tested_Allele", "BETA")], file = "bmi_snp_mom_model3.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



