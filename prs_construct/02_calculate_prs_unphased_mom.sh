# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 02_calculate_prs_unphased_mom
# Description: Calculate maternal Polygenic Risk Scores (PRS) for BMI 
#              using PLINK in the ALSPAC Mother-only cohort  
# ------------------------------------------------------------------

# ------------------------------
# Run PLINK to calculate PRS score mom (UKBB)
# ------------------------------
# Original PRS (model 1)
plink --bfile filtered_all.B3723.v3.MotherOnly \
      --score bmi_snp_mom_model1.txt 1 2 3 \
      --out prs_mom_model1

# # PRS without SNPs associated with T2D related traits (model 2)
plink --bfile filtered_all.B3723.v3.MotherOnly \
      --score bmi_snp_mom_model2.txt 1 2 3 \
      --out prs_mom_model2

# PRS without child BMI SNPs (model 3)
plink --bfile filtered_all.B3723.v3.MotherOnly \
      --score bmi_snp_mom_model3.txt 1 2 3 \
      --out prs_mom_model3

# PRS without pleiotropic SNPs (model 4)
plink --bfile filtered_all.B3723.v3.MotherOnly \
      --score bmi_snp_mom_model4.txt 1 2 3 \
      --out prs_mom_model4


