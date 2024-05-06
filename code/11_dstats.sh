### ============================================================================
### phylo2
### D-statistics
### ============================================================================

### Preparations

mkdir /gss/work/haex1482/phylo2/8_dstats

base=/gss/work/haex1482/phylo2/8_dstats && cd $base

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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-env/8.3
ml VCFtools/0.1.16-GCC-8.3.0
ml BCFtools/1.15.1-GCC-8.3.0

vcf=/gss/work/haex1482/phylo2/2_genotyping/out/9_filt

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
bgzip > phyps2e_m2n1.vcf.gz

tabix -p vcf phyps2e_m2n1.vcf.gz

bcftools +prune -m 0.5 -w 50kb -n 100 phyps2e_m2n1.vcf.gz -Oz -o phyps2e_m2n1l5.vcf.gz
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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-uniol-new-env libstdc++/7.1.0

gzip -cd phylo2e_m2n1l5.vcf.gz > phylo2e_m2n1l5.vcf

Dsuite Dtrios \
  phylo2e_m2n1l5.vcf \
  sets_phylo2e.tsv \
  -o phylo2e_m2n1l5.dtrios \
  -c

rm phylo2e_m2n1l5.vcf

# gzip -cd phyps2e_m2n1l5.vcf.gz > phyps2e_m2n1l5.vcf

# Dsuite Dquartets \
#   phyps2e_m2n1l5.vcf \
#   sets_phyps2e.tsv \
#   -o phyps2e_m2n1l5.dquart

# rename dquart_quartets dquarts *
# rm phyps2e_m2n1l5.vcf


### ============================================================================
### 3. D-stats analysis and plotting

# phylo2e Dtrios: 75 ingroup populations, 67525 trios

# see dstats_phylo2e.R

# m2n1l5 filter, Bonferroni correction:
# Dmin: x trios (y%) with significant deviation of D
# BBAA: x trios (y%) with significant deviation of D

## Plotting with Ruby (run locally from 5_Phylogeny/dstats)
# ruby ../code/plot_d.rb \
#   phylo2e_ld5n.dtrios_BBAA.txt \
#   order_alpha.txt \
#   0.02 \
#   heatmap_D_phylo2e_ld5n.dtrios_BBAA.svg

## Plotting in R: see R/dstats_phylo2e.R


### ============================================================================
### 4. Fbranch

mkdir /fs/dss/work/haex1482/phylo2/8_dstats/fbranch

base=/fs/dss/work/haex1482/phylo2/8_dstats/fbranch && cd $base

### ----------------------------------------------------------------------------
### Dataset

for i in tab tor atl eco casliz gemflk gumhon libhai nigboc maybel prohon
do
  grep "$i" /fs/dss/work/haex1482/phylo2/0_metadata/ids_phylo2e.txt |
    sort -R |
    head -n 6 |
    sort \
    >> ids_mono1.txt
done

# mono1
# 20478tabhon
# 20480tabhon
# 28393torboc
# Bocas16.3torboc
# Bocas16.4torboc
# s_tort_3torboc
# 52988atltam
# 52989atltam
# 52990atltam
# 54689atlliz
# 54761atlliz
# 54786atlliz
# 62549ecoarc
# 62556ecoarc
# 62576ecoarc
# 54649casliz
# 54650casliz
# PL17_142gemflk
# PL17_144gemflk
# PL17_145gemflk
# PL17_148gemflk
# PL17_153gemflk
# 20420gumhon
# 20426gumhon
# 20615gumhon
# 20617gumhon
# 20642gumhon
# 20643gumhon
# HypoHaiti1libhai
# HypoHaiti2libhai
# HypoHaiti3libhai
# 18418nigboc
# 18901nigboc
# 18904nigboc
# 18905nigboc
# 18906nigboc
# 18907nigboc
# PL17_119maybel
# PL17_120maybel
# PL17_121maybel
# PL17_122maybel
# PL17_123maybel
# PL17_126maybel
# 20650prohon
# 20845prohon
# 20846prohon

# mono2 adds:
# PL17_134uniflk
# PL17_135uniflk
# PL17_136uniflk
# PL17_140uniflk
# PL17_141uniflk
# PL17_143uniflk
# 62953puearc
# 62559floarc


### ----------------------------------------------------------------------------

#!/bin/bash

#SBATCH --job-name=remove
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00-02  # DD-HH
#SBATCH --output=sl_%j_remove.err
#SBATCH --error=sl_%j_remove.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml VCFtools/0.1.16-GCC-13.1.0
ml BCFtools/1.18-GCC-13.1.0

base=/fs/dss/work/haex1482/phylo2/8_dstats/fbranch

vcftools \
    --gzvcf $base/../../2_genotyping/out/9_filt/phylo2e_snpsfilt.vcf.gz \
    --keep ids_mono2.txt \
    --mac 2 \
    --max-missing 1 \
    --recode \
    --stdout |
grep -v -e '##contig' -e '##GATKCommandLine=' |
bgzip > mono2_m2n1.vcf.gz
#> 46 samples, 16037952 sites

tabix -p mono2_m2n1.vcf.gz

bcftools +prune -m 0.5 -w 50kb -n 100 mono2_m2n1.vcf.gz -Oz -o mono2_m2n1l5.vcf.gz
#> mono1: 46 samples, 915518 sites
#> mono2: 54 samples, 922531 sites

### ----------------------------------------------------------------------------

# Convert to Nexus format
gzip -d mono2_m2n1l5.vcf.gz > mono2_m2n1l5.vcf

ruby ~/apps/convert_vcf_to_nexus.rb mono2_m2n1l5.vcf mono2_m2n1l5.nex

# Obtain outgroup taxon numbers
grep '#CHROM' mono2_m2n1l5.vcf |
  sed 's/\s/\n/g' |
  tail -n +10 |
  nl \
  > sample_numbers2.txt

grep 'tor\|tab' sample_numbers.txt
#  9	20478tabhon
# 10	20480tabhon
# 18	28393torboc
# 30	Bocas16.3torboc
# 31	Bocas16.4torboc
# 46	s_tort_3torboc

grep 'tor\|tab' sample_numbers2.txt
#  9	20478tabhon
# 10	20480tabhon
# 18	28393torboc
# 32	Bocas16.3torboc
# 33	Bocas16.4torboc
# 54	s_tort_3torboc

# gzip mono2_m2n1l5.vcf
