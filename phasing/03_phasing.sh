# ------------------------------------------------------------------
# bash version 3.2.57
# File name: 03_phasing
# Description: phasing using phase_common function 
# ------------------------------------------------------------------

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-6
#SBATCH --output=log/03_phasing.%A_%a.out
#SBATCH --error=log/03_phasing.%A_%a.err

module load shapeit5
module load bcftools/1.22 bcftools

CHR="$((SLURM_ARRAY_TASK_ID))"

MAP=map/chr${CHR}.b37.gmap.gz
INP=input/filtered_all.B3723.v3.chr${CHR}.test.vcf.gz
OUT=output/filtered_all.B3723.v3.chr${CHR}.phased.bcf
LOG=log/filtered_all.B3723.v3.chr${CHR}.phased.log
PED=pedigree.txt

phase_common --input $INP --map $MAP --output $OUT --pedigree $PED --region ${CHR} --log $LOG --thread 8

bcftools index -f $OUT --threads 8