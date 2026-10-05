# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 03_analysis_with_puberty_phased
# Description: TSLS and TSRI with fixed covariates 
#              first stage: AgeDelivery_Mom, top 3 pcs 
#              second stage: sex, pubertal age (scaled)
# ------------------------------------------------------------------

library(boot)

# All outcomes and timepoints  
outcomes <- c("zfglucose_", "zfinsulin_", "zhomair_", "dysglycemie_", "InsulinResistance_")
time_points <- c("F7", "F15", "F18", "F24")
instruments <- c("SCORE_mom")

# Mapping of outcome variable prefixes to display labels
prefix_labels <- list(
  zfglucose_ = "Fasting Glucose",
  zfinsulin_ = "Fasting Insulin",
  zhomair_ = "HOMA-IR",
  dysglycemie_ = "Dysglycemia",
  InsulinResistance_ = "Insulin Resistance"
)

top3_pc_mom <- c("PC1_mom", "PC2_mom", "PC3_mom") # PCs 
first_covs <- c("AgeDelivery_Mom", top3_pc_mom) # covariates in first stage 
second_covs <- c("sex") # covariates in second stage 

# Data frame to store results 
results_table <- data.frame(
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

# Iterate through each outcome 
for (selected_outcome in outcomes) {
  
  # Retrieve the display label for the selected outcome
  outcome_label <- prefix_labels[[selected_outcome]]

  # Iterate through each time point  
  for (tp in time_points) {

    # Create visit-specific variables 
    second_covs_tp <- second_covs

    # Add puberty adjustment only at age 15
    if (tp == "F15") {
      second_covs_tp <- c(second_covs_tp, "puberty")
    }

    outcome_var <- paste0(selected_outcome, tp)
    age_label <- tp 
    
    # Retain participants with complete data for the exposure, outcome,
    # genetic instruments, and covariates required for the analysis.
    all_vars <- c("prepreg_bmi", outcome_var, instruments, second_covs_tp, first_covs)
    valid <- complete.cases(pheno_prs[, all_vars])
    pheno_valid <- pheno_prs[valid, ]
    
    # Perform TSLS or TSRI, depending whether the outcome is binary or continuous 
    if (selected_outcome == "dysglycemie_" | selected_outcome == "InsulinResistance_") {

      # Bootstrap 
      set.seed(42) # For reproducibility
      boot_tsri <- boot(
        data = pheno_valid,
        statistic = tsri,
        instruments = instruments,
        first_covs = first_covs,
        second_covs = second_covs_tp,
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
        instruments = instruments,
        first_covs = first_covs,
        second_covs = second_covs_tp,
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
    if (!is.na(n_cases)) { N_display <- paste0(N, " (", n_cases, ")") } else { N_display <- as.character(N) } # formatting N 
  
    # Append to results table
    results_table <- rbind(results_table, data.frame(
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

print(significant_results)
write.csv(results_table, file = paste0('model', model, "_with_puberty_phased.csv"), row.names = FALSE)
