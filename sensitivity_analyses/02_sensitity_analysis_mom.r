# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 02_sensitivity_analysis_mom
# Description: Sensitivity analysis #3 on the female-specific bmi-associated SNPs
# ------------------------------------------------------------------

# --------------------------
# Load Libraries
# --------------------------
library(TwoSampleMR) 
library(data.table)
library(usethis)
library(devtools)
library(biomaRt)
library(BSgenome)
library(dplyr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(MRPRESSO)
library(ggplot2)
library(dplyr)
library(patchwork)

# --------------------------
# Exposure data
# --------------------------
exposure <- read.table("exposure_data_mom_sensitivity_analysis_3.txt", header = TRUE)

# --------------------------
# Outcome data: 
# t2d_outcome = adult T2D outcome - Hg19  
# --------------------------
freadf <- function(df){ return(as.data.frame(fread(df))) }  # Helper function
t2d_outcome <- freadf("adult_t2d_gwas.txt")

t2d_outcome$chrposID <- sub("^chr", "", t2d_outcome$chrposID)

# Renaming columns (T2D)
t2d_outcome <- t2d_outcome %>%
    rename(
        SNP = chrposID,
        beta.outcome = 'Fixed-effects_beta',
        effect_allele.outcome = effect_allele,
        other_allele.outcome = other_allele,
        eaf.outcome = effect_allele_frequency,
        se.outcome = 'Fixed-effects_SE',
        pval.outcome = 'Fixed-effects_p-value'
    )

# Creating columns (T2D)
t2d_outcome$id.outcome <- 'T2D'
t2d_outcome$outcome <- 'T2D'

# --------------------------
# Harmonisation and Two-Sample MR
# --------------------------

# Harmonising data 
data <- harmonise_data(
  exposure_dat = exposure,
  outcome_dat = t2d_outcome
) 

write.table(data, "harmonised_mom_data.txt", row.names = FALSE)

# MR Egger and pleiotropy test 
results <- mr(data, method_list = c("mr_egger_regression"))
pleiotropy_test <- mr_pleiotropy_test(data)

knitr::kable(results)
knitr::kable(pleiotropy_test)

# LOO analysis 
res_loo <- mr_leaveoneout(data)

# Plot: too many SNPs so print half of the SNPs side by side 
n <- nrow(res_loo)

res1 <- res_loo[1:ceiling(n/2), ]
res2 <- res_loo[(ceiling(n/2) + 1):n, ]

p1 <- mr_leaveoneout_plot(res1)[[1]]
p2 <- mr_leaveoneout_plot(res2)[[1]]

combined_loo <- p1 | p2

# MR PRESSO
presso_results <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE,
          data = data)

# Heterogeneity test
heterogeneity_test <- mr_heterogeneity(data)


