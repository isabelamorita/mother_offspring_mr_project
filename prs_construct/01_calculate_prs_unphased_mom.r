# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_calculate_prs_unphased_mom
# Description: Construct a polygenic risk score (PRS) for BMI in females
#               and find the best set of parameters to find the PRS that
#               that explains more variance in the maternal BMI data. 
# ------------------------------------------------------------------

# --------------------------
# Load Libraries
# --------------------------
library(TwoSampleMR)
library(dplyr)
library(data.table)
library(genetics.binaRies)
library(rtracklayer)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(reshape2)

# --------------------------
# Load Phenotype Data
# --------------------------
phenotype <- readRDS('processed_pheno.rds')
phenotype$IID <- paste0(phenotype$gi_hrc_g0m_g1, "M")  # Add 'M' suffix for mothers

# Add PCA
pca <- read.table('PC.MotherOnly', header = TRUE)
colnames(pca)[3:20] <- paste0(colnames(pca)[3:20], "_mom")  # Append _mom to principal component names

# Merge phenotype data with PCA data using IID
phenotype <- merge(phenotype, pca[, 2:ncol(pca)], by = 'IID', all.x = TRUE)

# ---------------------
# Load GWAS summary statistics
# ---------------------
freadf <- function(df){ return(as.data.frame(fread(df))) } # Helper function: Load file as data.frame

# Load and preprocess GWAS exposure summary statistics (BMI in adult females)
exposure <- freadf('bmi.giant-ukbb.meta-analysis.females.txt')
exposure$SNP <- sub(':.*', '', exposure$SNP)  # Remove annotations after colon
names(exposure)[names(exposure)=='CHR'] <- "chr_name"
names(exposure)[names(exposure)=='POS'] <- "chrom_start"
names(exposure)[names(exposure)=='P'] <- "pval.exposure"
exposure_original <- exposure  # Save original version for reuse

# Define parameter thresholds for PRS generation
pval <- 5e-8
fstat <- 10
clump_r2 <- 0.001

# Initialize empty dataframe to store PRS regression results
results_table <- data.frame()


# ---------------------
# Construct PRS
# ---------------------
# Step 1: Filter exposure data based on p-value, remove rare variants and HLA 
exposure_filt <- exposure_original
exposure_filt <- exposure_filt %>% filter(!(chr_name == 6 & chrom_start >= 28510120 & chrom_start <= 33480577)) # Remove HLA
exposure_filt <- exposure_filt[nchar(exposure_filt$Tested_Allele)== 1 &
               nchar(exposure_filt$Other_Allele)== 1 &
                  exposure_filt$Tested_Allele %in% c("A" , "T", "C" , "G") &
          exposure_filt$Other_Allele %in% c("A" ,"T", "C","G"),] # Remove indels
exposure_filt <-subset(exposure_filt,Freq_Tested_Allele>0.01 & Freq_Tested_Allele<0.99) # Remove rare variants 
exposure_filt$p_value <- as.numeric(exposure_filt$pval.exposure)
exposure_filt <- exposure_filt[exposure_filt$p_value < pval, ] # p-value filtering

# Calculate R2 and F-statistic - Filter based on F-stat 
exposure_filt$R2 <- get_r_from_bsen(exposure_filt$BETA, exposure_filt$SE, exposure_filt$N)^2
exposure_filt$Fstat <- (exposure_filt$R2 / 1) / ((1 - exposure_filt$R2) / (exposure_filt$N - 2))
exposure_filt <- exposure_filt[exposure_filt$Fstat > fstat,]

# Step 2: LD clumping to remove correlated SNPs
clumped_all <- data.frame()
for (chr in 1:22) {
  chr_data <- exposure_filt[exposure_filt$chr_name == chr, ]
  bfile_path <- paste0('Chr', chr, '.EUR')
  
  tryCatch({
    clumped_chr <- clump_data(
      dat = chr_data,
      clump_kb = 10000,
      clump_r2 = clump_r2,
      clump_p1 = 1,
      clump_p2 = 1,
      pop = "EUR",
      bfile = bfile_path,
      plink_bin = (genetics.binaRies::get_plink_binary())
    )
    clumped_all <- rbind(clumped_all, clumped_chr)
  }, error = function(e) {
    cat("Clumping failed on chr", chr, ":", e$message, "\n")
  })
}

# Skip if not enough SNPs after clumping
if (nrow(clumped_all) < 10) {
  cat("Too few SNPs after clumping. Skipping...\n")
  next
}

# Create variant ID 
clumped_all$variant <- paste(clumped_all$chr_name, clumped_all$chrom_start, sep = ":")

# Step 3: Write PRS score file for PLINK
prs_file <- paste0("PRS_BMI_p", pval, "_F", fstat, "_r2", clump_r2, ".txt")
score_file <- clumped_all[, c("variant", "Tested_Allele", "BETA")]
write.table(score_file, prs_file, quote = FALSE, col.names = FALSE, row.names = FALSE)
num_snps <- nrow(score_file)

# Step 4: Run PLINK to compute PRSs
out_prefix <- gsub(".txt", "_output", prs_file)
plink_cmd <- paste(
  "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/plink/2.00-20231024-avx2/plink2",
  "--bfile filtered_all.B3723.v3.MotherOnly",
  "--score", prs_file, "1 2 3",
  "--out", out_prefix
)
system(plink_cmd)

# Step 5: Read in PLINK scores and merge with phenotype data
score_path <- paste0(out_prefix, ".sscore")

if (!file.exists(score_path)) next
prs_out <- read.delim(score_path, header = TRUE)
colnames(prs_out)[1] <- "FID"
colnames(prs_out)[5] <- "SCORE"
prs_out$SCORE <- scale(prs_out$SCORE)  # Standardize scores

pheno_prs <- merge(phenotype, prs_out[, c("FID", "IID", "SCORE")], by = "IID")

# Step 6: Linear regression of BMI on PRS
model <- lm(prepreg_bmi ~ SCORE, data = pheno_prs)
prs.coef <- summary(model)$coeff["SCORE", ]
prs.result <- data.frame(
  PRS_file = prs_file,
  P_threshold = pval,
  Fstat_threshold = fstat,
  clump_r2 = clump_r2,
  Num_SNPs = num_snps,
  R2 = summary(model)$r.squared,
  Beta = prs.coef[1],
  SE = prs.coef[2],
  P = prs.coef[4],
  F_stat = summary(model)$fstatistic[1]
)

# Save final results summary
write.csv(prs.result, "prs_result.csv", row.names = FALSE)
print(prs_file)