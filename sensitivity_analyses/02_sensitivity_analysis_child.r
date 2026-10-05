# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 02_sensitivity_analysis_child
# Description: Sensitivity analysis #3 on the child-specific bmi-associated SNPs
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

# --------------------------
# Exposure data
# --------------------------
exposure <- read.table("exposure_data_child_sensitivity_analysis_3.txt", header = TRUE, sep = "\t",)

# --------------------------
# Outcome data: 
# t2d_outcome = child T2D outcome - Hg19  
# --------------------------
freadf <- function(df){ return(as.data.frame(fread(df))) }  # Helper function
t2d_outcome <- freadf("child_t2d_gwas.txt")

# Renaming columns (T2D)
t2d_outcome <- t2d_outcome %>%
    rename(
        beta.outcome = BETA,
        info = SNP,
        effect_allele.outcome = EFFECT_ALLELE,
        other_allele.outcome = NON_EFFECT_ALLELE,
        eaf.outcome = EFFECT_ALLELE_FREQ,
        pval.outcome = P_VAL
    )

# Creating columns (T2D)
t2d_outcome$SNP <- paste0(t2d_outcome$Chr, ":", t2d_outcome$Position) # SNP 
t2d_outcome$t_val<-qt(t2d_outcome$pval.outcome/2,df=98) 
t2d_outcome$se.outcome<-t2d_outcome$beta.outcome/abs(t2d_outcome$t_val)
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

write.table(data, "harmonised_child_data.txt", row.names = FALSE) 

# MR Egger and pleiotropy test 
results <- mr(data, method_list = c("mr_egger_regression"))
pleiotropy_test <- mr_pleiotropy_test(data)

# LOO analysis 
res_loo <- mr_leaveoneout(data)
LOO_plot <- mr_leaveoneout_plot(res_loo)
LOO_plot[[1]]
ggsave(LOO_plot[[1]], file = "LOO_plot_downie.png", width = 7, height = 7)

# Filtering for only positive standard error values 
# Necessary for MR PRESSO
data <- data[data$se.exposure > 0 & data$se.outcome > 0, ]

# MR PRESSO
presso_results <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE,
          data = data)

# Heterogeneity test
heterogeneity_test <- mr_heterogeneity(data)
