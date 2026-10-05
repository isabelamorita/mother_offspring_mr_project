# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_sensitivity_analysis_child
# Description: Sensitivity analysis #1 on the child-specific bmi-associated SNPs
# ------------------------------------------------------------------

# --------------------------
# Load Libraries
# --------------------------
library(rtracklayer)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(reshape2)
library(data.table)
library(TwoSampleMR)
library(dplyr)
library(genetics.binaRies)


# Child BMI-associated SNP list (with columns chr_name, chr_position, effect_allele, other_allele and effect_weight)
child_bmi_snps <- read.table("child_bmi_snp_file.txt", comment.char = "#", header = TRUE)
child_bmi_snps$variant <- paste(child_bmi_snps$chr_name, child_bmi_snps$chr_position, sep=":")

# --------------------------
# Outcome GWAS 
# --------------------------
freadf <- function(df){return(as.data.frame(fread(df)))} # Helper function: Load file as data.frame

# Filter 1
fg_outcome <- freadf("FG_EUR_GCST90002232.tsv")
fi_outcome <- freadf("FI_EUR_GCST90002238.tsv")
hba1c_outcome <- freadf("HbA1c_EUR_GCST90002244.tsv")
homair_outcome <- freadf("HOMAIR_GCST90267567.tsv")
# homair_outcome is in GRCh38, the rest is in GRCh37

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

# Filtering for fasting glucose, fasting insulin and hba1c 
child_bmi_snps <- child_bmi_snps[!child_bmi_snps$variant %in% fg_outcome$variant, ]
child_bmi_snps <- child_bmi_snps[!child_bmi_snps$variant %in% fi_outcome$variant, ]
child_bmi_snps <- child_bmi_snps[!child_bmi_snps$variant %in% hba1c_outcome$variant, ]

# ---------------------
# Liftover of PRS variant coordinates from hg19 to hg38
# ---------------------

chain <- import.chain("hg19ToHg38.over.chain")

# Add "chr" prefix to chromosomes 
child_bmi_snps$chr_name <- ifelse(grepl("^chr",child_bmi_snps$chr_name),
                        child_bmi_snps$chr_name,
                        paste0("chr",child_bmi_snps$chr_name))

# Ensure base pair positions are numeric
child_bmi_snps$chr_position <- as.numeric(child_bmi_snps$chr_position)

# Create a GRanges object for liftover
gr <- GRanges(
  seqnames = child_bmi_snps$chr_name,
  ranges = IRanges(start = child_bmi_snps$chr_position, end = child_bmi_snps$chr_position),
  chromosome = child_bmi_snps$chr_name, 
  base_pair_location = child_bmi_snps$chr_position, 
  effect_allele = child_bmi_snps$effect_allele, 
  beta = child_bmi_snps$effect_weight, 
)
print(head(gr))

# Perform the liftover
lifted_hg38 <- liftOver(gr, chain) # Perform conversion
lifted_gr_hg38 <- unlist(lifted_hg38) # Flatten result into a GRanges object
lifted_df_hg38 <- as.data.frame(lifted_gr_hg38)

# Remove 'chr' prefix for compatibility with downstream datasets
lifted_df_hg38$chromosome <- sub("^chr", "", lifted_df_hg38$chromosome) # Removing the 'chr' prefix from chromosome number

# Create a new variant ID in the format CHR:POS
lifted_df_hg38$variant_id <- paste0(lifted_df_hg38$chromosome, ":", lifted_df_hg38$start) 

# Filtering for homa-ir snps 
homair_outcome$variant_id <- sub("^([^_]+)_([^_]+)_.*$", "\\1:\\2", homair_outcome$variant)
lifted_df_hg38 <- lifted_df_hg38[!lifted_df_hg38$variant_id %in% homair_outcome$variant_id, ] # no homa-ir snps removed 

# FILTER 1
write.table(child_bmi_snps, file = "bmi_snp_child_model2.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
