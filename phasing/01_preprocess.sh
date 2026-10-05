# ------------------------------------------------------------------
# bash version 3.2.57
# File name: 01_preprocess
# Description: # Convert PLINK files to VCF (because SHAPEIT5 takes that as an input)
# ------------------------------------------------------------------

module load plink
module load bcftools/1.22 bcftools
module load htslib/1.22.1 bcftools

for CHR in $(seq 3 22); do
  plink --bfile filtered_all.B3723.v3 \
        --chr $CHR \
        --recode vcf-iid \
        --out /input/filtered_all.B3723.v3.chr${CHR}
  bgzip -f /input/filtered_all.B3723.v3.chr${CHR}.vcf
  bcftools index -t filtered_all.B3723.v3.chr${CHR}.vcf.gz
done

# Download map files (GRCh37) and extract them 
# link from github: https://github.com/odelaneau/shapeit5/blob/main/resources/maps/b37/genetic_maps.b37.tar.gz
# extract: tar -xvzf genetic_maps.b37.tar.gz


