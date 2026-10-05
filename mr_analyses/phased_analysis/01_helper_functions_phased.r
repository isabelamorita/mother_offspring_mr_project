# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_helper_functions_phased 
# Description: Functions to (1) load and preprocess the data and (2) perform MR analysis 
# ------------------------------------------------------------------

library(haven)

# Function to load phenotype data with PRS data 
load_pheno_with_prs <- function(model) {
  mom_files <- c(
    'prs_phased_genotype_mom_model1.txt',
    'prs_phased_genotype_mom_model2.txt',
    'prs_phased_genotype_mom_model3.txt',
    'prs_phased_genotype_mom_model4.txt'
  )
  
  # Pick files 
  mom_file <- mom_files$female[[model]]

  # Processed phenotype file
  processed_pheno <- readRDS('processed_pheno.rds')
  
  # Mom
  PRS_mom <- read.table(mom_file, header=T)
  PRS_mom$SCORE<-scale(PRS_mom$SCORE) 
  PRS_mom$IID_mom <- PRS_mom$IID
  PRS_mom$SCORE_mom <- PRS_mom$SCORE
  pheno_prs <- merge(processed_pheno, PRS_mom[,c('FID','IID_mom','SCORE_mom')], by="IID_mom")

  return(pheno_prs)
}

# TSLS: MR analysis for continuous outcomes 
tsls <- function(data, instruments, first_covs, second_covs, outcome_var){

    # First stage
    first_formula <- as.formula(paste("prepreg_bmi ~", paste(c(instruments, first_covs), collapse = " + ")))
    first_stage <- lm(first_formula, data = data)
    data$bmi_fitted <- fitted(first_stage)
    
    # Second stage
    second_formula <- as.formula(
    paste(outcome_var, "~ bmi_fitted",
        if (length(second_covs) > 0)
        paste("+", paste(second_covs, collapse = " + "))
        else "")
    )
    
    second_stage <- lm(second_formula, data = data)
    
    return(second_stage)
}

# TSRI: MR analysis for binary outcomes
tsri <- function(data, indices, instruments, first_covs, second_covs, outcome_var){
    d <- data[indices, ]
    
    first_formula <- as.formula(paste("prepreg_bmi ~", paste(c(instruments, first_covs), collapse = " + ")))
    first_stage <- lm(first_formula, data = d)
    d$bmi_fitted <- fitted(first_stage)
    
    d$res <- residuals(first_stage)
          
    # Second stage with residual inclusion, logistic regression
    second_formula<- paste(outcome_var, "~ prepreg_bmi + res", if(length(second_covs) > 0) paste("+", paste(second_covs, collapse = " + ")) else "")
    second_stage <- glm(as.formula(second_formula), data = d, family = binomial(link = "logit"))
    
    return(coef(second_stage)["prepreg_bmi"])

}
