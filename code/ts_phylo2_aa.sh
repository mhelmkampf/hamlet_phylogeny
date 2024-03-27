### ============================================================================
### phylo2
### Inferring genomic history based on ancestral recombination graph (tskit)
### ============================================================================

### Preparations

base=$WORK/phylo2

mkdir $base/5_arg

cd $base/5_arg

mkdir 1_counts
mkdir 2_sfs
mkdir 3_vcf
mkdir 4_samples
mkdir 5_tsd
mkdir 6_pi
mkdir 7_tmrca
mkdir log

## -----------------------------------------------------------------------------
## Create outgroup id files

grep -E 'tig' $base/0_metadata/ids_phylo2e.txt > 1_counts/ids_tig.txt
grep -E 'tor|tab' $base/0_metadata/ids_phylo2e.txt > 1_counts/ids_tortab.txt

cat > 2_sfs/seedfile.txt <<EOF
96577627
EOF

cat > 2_sfs/config.txt <<EOF
n_outgroup 2
model 1
nrandom 1
EOF

cat > 3_vcf/HypPue1.1_LG.sizes <<EOF
LG01 25890264
LG02 31702000
LG03 23122023
LG04 25000711
LG05 20912025
LG06 24210077
LG07 24100755
LG08 27000498
LG09 27186781
LG10 21876389
LG11 21320898
LG12 25506241
LG13 26810734
LG14 22338278
LG15 21633613
LG16 24734982
LG17 28429713
LG18 19773282
LG19 27051319
LG20 24942031
LG21 20142184
LG22 13699534
LG23 17761902
LG24 14503443
LG_M 17159
EOF

## -----------------------------------------------------------------------------
## Notes on regions

# LG=LG12
# LENGTH=25506241
# START=20322000
# END=20358000
# PREFIX=${lg}-casz1c
# casz1:  LG12, 20135000-20300000 (200 kb)
# casz1b: LG12, 20220000-20250000 ( 30 kb, fst peak)
# casz1c: LG12, 20322000-20358000 ( 36 kb, gene cds)


### ============================================================================
### 0. Phasing (Shapeit)

# see $base/3_phasing: phylo2e_phased.vcf.gz
# (no mac, missingness or distance filter)


### ============================================================================
### 1. Infer ancestral alleles

#!/bin/bash

#SBATCH --job-name=aa
#SBATCH --partition=rosa_express.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --time=00-02   # DD-HH
#SBATCH --output=/dev/null
#SBATCH --error=log/sl_%A_%a_aa.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml VCFtools/0.1.16-GCC-13.1.0
ml GSL/2.7-GCCcore-13.1.0

base=$WORK/phylo2

lg=($(seq -w 1 24))


### Count alleles

vcftools \
    --gzvcf $base/4_phasing/4_merge/phylo2e_phased.vcf.gz \
    --remove 1_counts/ids_tortab.txt \
    --remove 1_counts/ids_tig.txt \
    --chr LG${lg[((SLURM_ARRAY_TASK_ID-1))]} \
    --counts \
    --stdout | gzip > 1_counts/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_ingroup.counts.gz

vcftools \
    --gzvcf $base/4_phasing/4_merge/phylo2e_phased.vcf.gz \
    --keep 1_counts/ids_tortab.txt \
    --chr LG${lg[((SLURM_ARRAY_TASK_ID-1))]} \
    --counts \
    --stdout | gzip > 1_counts/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_tortab.counts.gz

vcftools \
    --gzvcf $base/4_phasing/4_merge/phylo2e_phased.vcf.gz \
    --keep 1_counts/ids_tig.txt \
    --chr LG${lg[((SLURM_ARRAY_TASK_ID-1))]} \
    --counts \
    --stdout | gzip > 1_counts/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_tig.counts.gz


### Reconstruct AA (Est-sfs)

python ~/apps/est-sfs_format.py \
    -i 1_counts/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_ingroup.counts.gz \
    -o 1_counts/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_tortab.counts.gz \
     1_counts/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_tig.counts.gz

rename LG${lg[((SLURM_ARRAY_TASK_ID-1))]}. LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_ *
mv LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_* 2_sfs
rm 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_config.file

est-sfs \
    2_sfs/config.txt \
    2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_est.infile \
    2_sfs/seedfile.txt \
    2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_sfs.txt \
    2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_pvalues.txt


### Add AA to VCF

vcftools \
    --gzvcf $base/4_phasing/4_merge/phylo2e_phased.vcf.gz \
    --chr LG${lg[((SLURM_ARRAY_TASK_ID-1))]} \
    --recode \
    --stdout > 3_vcf/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}.vcf

sed '1,8d' 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_pvalues.txt |
awk '{ if ($4 > ($5 + $6 + $7)) print "A" ;
    else if ($5 > ($4 + $6 + $7)) print "C" ;
    else if ($6 > ($4 + $5 + $7)) print "G" ;
    else if ($7 > ($4 + $5 + $6)) print "T" }' \
    > 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.txt

paste --delimiters='\t' 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_pos.txt 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.txt \
    > 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.tsv && rm 2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.txt

vcf-info-annotator \
    3_vcf/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}.vcf \
    2_sfs/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.tsv \
    AA \
    -f Character \
    -d "Ancestral Allele" \
    -o 3_vcf/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.vcf

l=$(grep "LG${lg[((SLURM_ARRAY_TASK_ID-1))]}" 3_vcf/HypPue1.1_LG.sizes | cut -f2 -d " ")

awk -v lg=LG${lg[((SLURM_ARRAY_TASK_ID-1))]} -v l=$l 'NR == 5 { print "##contig=<ID="lg",length="l">" }1' 3_vcf/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.vcf \
    > 3_vcf/phylo2e_LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aal.vcf && \
    rm 3_vcf/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}.vcf 3_vcf/LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aa.vcf


### ============================================================================
### 2. Prepare final VCF files (no missing data)

#!/bin/bash

#SBATCH --job-name=filter
#SBATCH --partition=rosa_express.p
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --time=00-02  # DD-HH
#SBATCH --output=/dev/null
#SBATCH --error=log/sl_%A_%a_filter.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

base=$WORK/phylo2/5_arg

ml VCFtools/0.1.16-GCC-13.1.0

# lg=($(seq -w 1 24))

#   vcftools \
#     --vcf 3_vcf/phylo2e_LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_aal.vcf \
#     --max-missing 1 \
#     --recode \
#     --stdout > 3_vcf/phylo2e_LG${lg[((SLURM_ARRAY_TASK_ID-1))]}_n1.vcf

# done


### ----------------------------------------------------------------------------
### Optional: remove outgroup

for lg in LG01 LG02 LG03 LG04 LG08 LG12
do

  vcftools \
    --vcf $base/3_vcf/phylo2e_${lg}_aal.vcf \
    --remove $base/3_vcf/outgroup.ids \
    --max-missing 1 \
    --recode \
    --stdout > $base/3_vcf/phylo2e-S_${lg}_n1.vcf

done
