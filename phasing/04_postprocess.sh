# ------------------------------------------------------------------
# bash version 3.2.57
# File name: 04_postprocess 
# Description: Extract only BMI-associated SNPs from all bcf files 
# ------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=log/04_postprocessing.%j.out
#SBATCH --error=log/04_postprocessing.%j.err

module load bcftools/1.22 bcftools

# Step 1: Merge all chromosomes into one big BCF and index it 
bcftools concat -Ob -o merged_phased.bcf output/filtered_all.B3723.v3.chr{1..22}.phased.bcf

bcftools index merged_phased.bcf

# Step 2: Extract only BMI-associated SNPs
bcftools view -T mom_bmi_snps.txt \
    merged_phased.bcf -Ob -o bmi_all.bcf

# Step 3: Index the final file
bcftools index bmi_all.bcf