# ------------------------------------------------------------------
# R version - 4.4.0
# File name: 01_calculate_prs_child
# Description: Calculate maternal Polygenic Risk Scores (PRS) for BMI 
#              using PLINK in the ALSPAC Child-only cohort 
# ------------------------------------------------------------------

# ------------------------------
# Run PLINK to calculate PRS score child  
# ------------------------------
# Original PRS (model 1)
plink --bfile filtered_all.B3723.v3.ChildOnly \
      --score bmi_snp_child_model1.txt 1 2 3 \
      --out prs_child_model1

# PRS without SNPs associated with T2D related traits (model 2 and 3)
plink --bfile filtered_all.B3723.v3.ChildOnly \
      --score bmi_snp_child_model2.txt 1 2 3 \
      --out prs_child_model2

# PRS without pleiotropic SNPs (model 4)
plink --bfile filtered_all.B3723.v3.ChildOnly \
      --score bmi_snp_child_model4.txt 1 2 3 \
      --out prs_child_model4

