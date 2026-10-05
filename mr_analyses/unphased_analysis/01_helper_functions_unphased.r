# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_helper_functions_unphased
# Description: Functions to (1) load and preprocess the data and (2) perform MR analysis 
# ------------------------------------------------------------------

# Function to load phenotype data with PRS and PCs data 
load_pheno_with_prs <- function(model) {
  mom_files <- c(
    'prs_mom_model1.profile',
    'prs_mom_model2.profile',
    'prs_mom_model3.profile',
    'prs_mom_model4.profile'
  )

  child_files <- c(
    'prs_child_model1.profile',
    'prs_child_model2.profile',
    'prs_child_model3.profile',
    'prs_child_model4.profile'
  )
  
  # Pick files 
  mom_file <- mom_files$female[model] 
  child_file <- child_files[model]
  
  # Read phenotype data 
  processed_pheno <- readRDS('/home/imorita/links/scratch/data/BMI/processed_pheno.rds')
  
  # Mom
  PRS_mom <- read.table(mom_file, header=T)
  PRS_mom$SCORE<-scale(PRS_mom$SCORE) 
  PRS_mom$IID_mom <- PRS_mom$IID
  PRS_mom$SCORE_mom <- PRS_mom$SCORE
  pheno_prs_temp <- merge(processed_pheno, PRS_mom[,c('FID','IID_mom','SCORE_mom')], by="IID_mom")

  # Child
  PRS_child <- read.table(child_file, header=T)
  PRS_child$SCORE<-scale(PRS_child$SCORE) 
  PRS_child$IID_child <- PRS_child$IID
  PRS_child$SCORE_child <- PRS_child$SCORE
  pheno_prs <- merge(pheno_prs_temp, PRS_child[,c('FID','IID_child','SCORE_child')], by="IID_child")

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

