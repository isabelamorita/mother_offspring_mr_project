# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 02_analysis_fixed_cov_unphased
# Description: TSLS and TSRI with fixed covariates 
#              first stage: AgeDelivery_Mom, top 3 pcs 
#              second stage: sex
#              sensitivity analysis: adjusting SCORE_child (and the child top 3 PCs)
# ------------------------------------------------------------------

library(boot)

# All outcomes and timepoints  
outcomes <- c("zfglucose_", "zfinsulin_", "zhomair_", "dysglycemie_", "InsulinResistance_")
time_points <- c("F7", "F15", "F18", "F24")
instruments <- c("SCORE_mom")

prefix_labels <- list(
  zfglucose_ = "Fasting Glucose",
  zfinsulin_ = "Fasting Insulin",
  zhomair_ = "HOMA-IR",
  dysglycemie_ = "Dysglycemia",
  InsulinResistance_ = "Insulin Resistance"
)

top3_pc_mom <- c("PC1_mom", "PC2_mom", "PC3_mom") # Mother PCs 
top3_pc_child <- c("PC1_child", "PC2_child", "PC3_child") # Child PCs 
second_covs <- c("sex") # covariates in second stage 
first_covs <- c("AgeDelivery_Mom", top3_pc_mom)  # covariates in first stage 

# Data frame to store results 
results_table <- data.frame(
  PRS_type = character(),
  OffspringBMI_adjusted = character(),
  Outcome = character(),
  Age = character(),
  pval = numeric(),
  beta = numeric(),
  OR = numeric(),
  SE = numeric(),
  ci_combined = character(),
  N = character(),
  stringsAsFactors = FALSE
)

  
model <- 4 # model number 
pheno_prs <- load_pheno_with_prs(model)

for (bmi_adj in c(FALSE, TRUE)) {
  
  # Temporary variables 
  instruments_tmp <- instruments
  first_covs_tmp <- first_covs
  second_covs_tmp <- second_covs

  # When adjusting for offspring BMI PRS, add it as a covariate 
  if (bmi_adj) {
    instruments_tmp <- c(instruments, "SCORE_child")
    second_covs_tmp <- c("SCORE_child", second_covs, top3_pc_child) 
    first_covs_tmp <- c(first_covs, top3_pc_child)
  }

  # Iterate through each outcome 
  for (selected_outcome in outcomes) {
    
    # Retrieve the display label for the selected outcome
    outcome_label <- prefix_labels[[selected_outcome]]

    # Iterate through each time point 
    for (tp in time_points) {
      
      outcome_var <- paste0(selected_outcome , tp)
      age_label <- tp 

      # Retain participants with complete data for the exposure, outcome,
      # genetic instruments, and covariates required for the analysis.
      all_vars <- c("prepreg_bmi", outcome_var, instruments_tmp, second_covs_tmp, first_covs_tmp)
      valid <- complete.cases(pheno_prs[, all_vars])
      pheno_valid <- pheno_prs[valid, ]

      # Perform TSLS or TSRI, depending whether the outcome is binary or continuous 
      if (selected_outcome == "dysglycemie_" | selected_outcome == "InsulinResistance_") {
        # Bootstrap 
        set.seed(42) # For reproducibility
        boot_tsri <- boot(
          data = pheno_valid,
          statistic = tsri,
          instruments = instruments_tmp,
          first_covs = first_covs_tmp,
          second_covs = second_covs_tmp,
          outcome_var = outcome_var,
          R = 1000
        )

        # Coefficients 
        beta_bmi <- boot_tsri$t0 # original estimate 
        se_bmi <- sd(boot_tsri$t) # bootstrap standard error
        ci_bmi <- quantile(boot_tsri$t, c(0.025, 0.975))
        R <- length(boot_tsri$t) # number of bootstrap replicates 
        z_bmi <- beta_bmi / se_bmi # necessary to calculate p-value 
        pval  <- 2 * (1 - pnorm(abs(z_bmi))) # p-value 
        or_bmi <- exp(beta_bmi) # Odds ratio 
        ci_combined <- paste0("[", sprintf("%.2f", exp(ci_bmi[1])), ", ", sprintf("%.2f", exp(ci_bmi[2])), "]" ) # confidence interval 

        # number of cases 
        counts <- table(pheno_valid[[outcome_var]])
        n_cases <- counts["1"]

      } else {
        # Continuous outcomes 
        second_stage <- tsls(data = pheno_valid,
          instruments = instruments_tmp,
          first_covs = first_covs_tmp,
          second_covs = second_covs_tmp,
          outcome_var = outcome_var
        )

        # coefficients
        summary_model <- summary(second_stage)
        beta_bmi <- summary_model$coefficients["bmi_fitted", 1]
        se_bmi <- summary_model$coefficients["bmi_fitted", 2]
        pval <- summary_model$coefficients["bmi_fitted", 4]
        ci_bmi <- confint(second_stage)["bmi_fitted", ]
        or_bmi <- NA
        ci_combined <- paste0( "[", sprintf("%.2f", ci_bmi[1]), ", ", sprintf("%.2f", ci_bmi[2]), "]")

        n_cases <- NA 
      }

      N <- nrow(pheno_valid)  # total sample size
      if (!is.na(n_cases)) { N_display <- paste0(N, " (", n_cases, ")") } else { N_display <- as.character(N) } # formatting N to include n_cases if outcome is binary 
  
      # Append to results table
      results_table <- rbind(results_table, data.frame(
        PRS_type = prs_type,
        OffspringBMI_adjusted = ifelse(bmi_adj, "Yes", "No"),
        Outcome = outcome_label,
        Age = age_label,
        pval = pval,
        beta = beta_bmi,
        OR = or_bmi,
        SE = se_bmi,
        ci_combined = ci_combined,
        N = N_display,
        stringsAsFactors = FALSE
      ))
    }
  }
}

print(results_table)
write.csv(results_table, file = paste0('model', model, "_fixed_cov_unphased.csv"), row.names = FALSE)

