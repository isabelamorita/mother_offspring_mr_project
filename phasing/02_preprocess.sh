# ------------------------------------------------------------------
# bash version 3.2.57
# File name: 02_preprocess
# Description: Add AC and AN files and index files 
# ------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --array=1-21
#SBATCH --output=log/02_preprocess.%A_%a.out
#SBATCH --error=log/02_preprocess.%A_%a.err

module load bcftools/1.22 bcftools

CHR="$((SLURM_ARRAY_TASK_ID))"

# Add AC and AN files 
bcftools +fill-tags input/filtered_all.B3723.v3.chr${CHR}.vcf.gz -Oz -o input/filtered_all.B3723.v3.chr${CHR}.test.vcf.gz -- -t AN,AC

# Indexing files (necessary for shapeit5)
bcftools index input/filtered_all.B3723.v3.chr${CHR}.test.vcf.gz
