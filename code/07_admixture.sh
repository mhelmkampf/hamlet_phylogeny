### ================================================================================
### Radiation with reproductive isolation in the near-absence of phylogenetic signal
### 07. Admixture analysis (ADMIXTURE)
### By Martin Helmkampf, last edited 2025-06-16
### ================================================================================

### Preparations

base=$WORK/phylo2

mkdir $base/6_admix

cd $base/6_admix

mkdir admcv
mkdir bed


### ============================================================================
### Calculate admixture proportions (admcv / Admixture)

#!/bin/bash

#SBATCH --job-name=admcv
#SBATCH --partition=rosa_express.p
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00-02   # DD-HH
#SBATCH --output=sl_%j_admcv.out
#SBATCH --error=sl_%j_admcv.err

base=$WORK/phylo2

gzip -cd $base/2_genotyping/out/9_filt/phyps-snp_LDfilt.vcf.gz |
    sed 's/LG//g' \
    > phyps-snp_LDfilt.ed.vcf

# Convert to BED format
plink --vcf phyps-snp_LDfilt.ed.vcf \
    --make-bed \
    --allow-extra-chr \
    --out bed/phyps-snp_LDfilt

# Run ADMIXTURE
for k in {1..12}
do
    admixture \
    --cv -j24 \
    bed/phyps-snp_LDfilt.bed $k > admcv/phyps-snp_LDfilt_k${k}.out
done

mv *.P admcv
mv *.Q admcv

# Print CV error
for k in {1..12}
do
    grep 'CV' admcv/phyps-snp_LDfilt_k${k}.out \
    >> CV_phyps-snp_LDfilt.out
done

# Add sample ids to proportions
for k in {1..12}
do
    paste -d " " $base/0_metadata/ids_phyps2e.txt admcv/phyps-snp_LDfilt.${k}.Q |
        sed 's/ $//g' \
        > AdmcvProp_phyps-snp_LDfilt_k${k}.tsv
done

# Plotted with R/admcv_phyps2.R

# CV error (K=1): 0.59175
# CV error (K=2): 0.58147
# CV error (K=3): 0.57689
# CV error (K=4): 0.57688
# CV error (K=5): 0.57609 <<<
# CV error (K=6): 0.57783
# CV error (K=7): 0.58136
# CV error (K=8): 0.58323
# CV error (K=9): 0.58668
# CV error (K=10): 0.59434
# CV error (K=11): 0.59585
# CV error (K=12): 0.61051
