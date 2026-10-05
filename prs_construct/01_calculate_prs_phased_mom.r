# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_calculate_prs_phased_mom 
# Description: Calculate the maternal PRS based on the non-transmitted allele
# ------------------------------------------------------------------

# Input: phased data vcf 
# Output: Calculated PRS for the mothers (non-transmitted alleles)

# Reading vcf file 
vcf <- read.table("bmi_all.vcf", comment.char = "#", header = FALSE, sep = "\t")
# rows: 439 SNPs (BMI-associated SNPs)
# columns: SNP ID, reference and alternative allele, all samples 

# Read the header line (starts with a single #)
header <- strsplit(grep("^#CHROM", readLines("bmi_all.vcf"), value = TRUE), "\t")[[1]]
colnames(vcf) <- header

# Identify all child and mom columns
all_child_cols <- grep("_f2A$", header)
all_mom_cols   <- grep("_f2M$", header) 

# Make sure they match by numeric ID
all_child_ids <- sub("A$", "", header[all_child_cols])
all_mom_ids <- sub("M$", "", header[all_mom_cols])

# Keep only pairs with both mom and child
common_ids <- intersect(all_child_ids, all_mom_ids) # 6089

# Filter columns and keep cleaned IDs
child_cols_keep <- all_child_cols[all_child_ids %in% common_ids] # 6089
mom_cols_keep <- all_mom_cols[all_mom_ids %in% common_ids] # 6089
child_ids_keep <- all_child_ids[all_child_ids %in% common_ids]
mom_ids_keep <- all_mom_ids[all_mom_ids %in% common_ids]

n_moms <- length(mom_cols_keep) 
n_snps <- nrow(vcf) 

# BMI-associated SNPs
bmi_snps <- read.table("bmi_snp_mom_model1.txt")
colnames(bmi_snps) <- c("SNP", "EffectAllele", "Beta")

# Vector to place all the scores of the PRS 
maternal_prs <- numeric(n_moms)

# Looping through each mother / child pair
for(pair in 1:n_moms) {
  
  child_col <- child_cols_keep[pair]
  mom_col <- mom_cols_keep[pair]

  # temporary vector for this mom's non-transmitted alleles
  non_transmitted_allele <- character(n_snps) 

  # Looping through each SNP
  for(i in 1:n_snps) {
    # get REF and ALT
    ref <- vcf[i, 4] # REF column
    alt <- vcf[i, 5] # ALT column
    
    # get phased genotypes
    child_gt <- strsplit(vcf[i, child_col], ":")[[1]][1]  # "0|1"
    mom_gt   <- strsplit(vcf[i, mom_col], ":")[[1]][1]    # "0|1"
    
    # transmitted allele (maternal) is right of "|"
    transmitted <- strsplit(child_gt, "\\|")[[1]][2]
    
    # mom alleles
    mom_alleles <- strsplit(mom_gt, "\\|")[[1]]
    
    if (mom_alleles[1] == transmitted) {
      non_transmitted <- mom_alleles[2]
    } else if (mom_alleles[2] == transmitted) {
      non_transmitted <- mom_alleles[1]
    } else {
      non_transmitted <- transmitted # both differ, possible phasing issue?
    }

    # convert to base (REF/ALT)
    if(non_transmitted == "0") {
      non_transmitted_allele[i] <- ref
    } else if(non_transmitted == "1") {
      non_transmitted_allele[i] <- alt
    } else {
      non_transmitted_allele[i] <- NA  # in case of missing or unexpected genotype
    }
  }

  # Match SNPs with the ones in 'bmi_snps'
  match_idx <- match(vcf$ID, bmi_snps$SNP)  # NA if SNP not in bmi_snps
  valid <- !is.na(match_idx) # only SNPs present in both

  # Extract vectors for PRS calculation
  allele_to_use <- non_transmitted_allele[valid]
  effect_allele <- bmi_snps$EffectAllele[match_idx[valid]]
  beta_to_use <- bmi_snps$Beta[match_idx[valid]]

  # Calculate PRS for this mom 
  prs <- sum((allele_to_use == effect_allele) * beta_to_use, na.rm = TRUE)
  maternal_prs[pair] <- prs

}

# Create PRS table
prs_df <- data.frame(
  FID = paste0(mom_ids_keep, "M"),
  IID = paste0(mom_ids_keep, "M"),
  SCORE = maternal_prs
)

# nt = non-trasmitted
write.table(prs_df, "prs_phased_genotype_mom_model4.txt", quote = FALSE, sep = "\t", row.names = FALSE)














