### ================================================================================
### Radiation with reproductive isolation in the near-absence of phylogenetic signal
### 11. D-Statistics (Dsuite)
### By Martin Helmkampf, last edited 2025-06-16
### ================================================================================

### Preparations

mkdir $WORK/phylo2/8_dstats

base=$WORK/phylo2/8_dstats && cd $base

# Create set files for Dsuite
sed -E 's/(.+)([a-z]{6})/\1\2\t\2/g' ../0_metadata/ids_phylo2e.txt |
  sed -E 's/\ttab.*|\ttig.*|\ttor.*/\tOutgroup/g' > sets_phylo2e.tsv

sed -E 's/(.+)([a-z]{6})/\1\2\t\2/g' ../0_metadata/ids_phyps2e.txt > sets_phyps2e.tsv


### ============================================================================
### 1. Prune VCF by linkage disequilibrium

#!/bin/bash

#SBATCH --job-name=prune
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --time=02-00
#SBATCH --output=/dev/null
#SBATCH --error=sl_%j_prune.err

ml hpc-env/8.3
ml VCFtools/0.1.16-GCC-8.3.0
ml BCFtools/1.15.1-GCC-8.3.0

vcf=$WORK/phylo2/2_genotyping/out/9_filt

vcftools \
    --gzvcf $vcf/phylo2e_snpsfilt.vcf.gz \
    --mac 2 \
    --max-missing 1 \
    --recode \
    --stdout |
grep -v -e '##contig' -e '##GATKCommandLine=' |
bgzip > phylo2e_m2n1.vcf.gz

tabix -p vcf phylo2e_m2n1.vcf.gz

bcftools +prune -m 0.5 -w 50kb -n 100 phylo2e_m2n1.vcf.gz -Oz -o phylo2e_m2n1l5.vcf.gz
#> 936095 sites (sites with r2 > 0.5 in 50 kb windows discarded)

vcftools \
    --gzvcf $vcf/phyps2e_snpsfilt.vcf.gz \
    --mac 2 \
    --max-missing 1 \
    --recode \
    --stdout | 
grep -v -e '##contig' -e '##GATKCommandLine=' |
bgzip > phyps-snpn1.vcf.gz

tabix -p vcf phyps-snpn1.vcf.gz

bcftools +prune -m 0.5 -w 50kb -n 100 phyps-snpn1.vcf.gz -Oz -o phyps-snp_LDfilt.vcf.gz
#> 935941 sites

# Move VCFs to 2_genotyping/out/9_filt
# mv *.vcf.gz* /fs/dss/work/haex1482/phylo2/2_genotyping/out/9_filt


### ============================================================================
### 2. Calculate D-stats (Dsuite)

#!/bin/bash

#SBATCH --job-name=dsuite
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --time=00-10
#SBATCH --output=/dev/null
#SBATCH --error=sl_%j_dsuite.err

ml hpc-uniol-new-env libstdc++/7.1.0

gzip -cd phylo2e_m2n1l5.vcf.gz > phylo2e_m2n1l5.vcf

Dsuite Dtrios \
  phylo2e_m2n1l5.vcf \
  sets_phylo2e.tsv \
  -o phylo2e_m2n1l5.dtrios \
  -c

rm phylo2e_m2n1l5.vcf
