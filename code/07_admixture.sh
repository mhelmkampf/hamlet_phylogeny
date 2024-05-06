### ============================================================================
### phylo2
### Estimate admixture proportions (for each k, per sample)
### ============================================================================

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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

base=$WORK/phylo2

gzip -cd $base/2_genotyping/out/9_filt/phyps2e_m2n1l5.vcf.gz |
    sed 's/LG//g' \
    > phyps2e_m2n1l5.ed.vcf

# Convert to BED format
plink --vcf phyps2e_m2n1l5.ed.vcf \
    --make-bed \
    --allow-extra-chr \
    --out bed/phyps2e_m2n1l5

# Run ADMIXTURE
for k in {1..12}
do
    admixture \
    --cv -j24 \
    bed/phyps2e_m2n1l5.bed $k > admcv/phyps2e_m2n1l5_k${k}.out
done

mv *.P admcv
mv *.Q admcv

# Print CV error
for k in {1..12}
do
    grep 'CV' admcv/phyps2e_m2n1l5_k${k}.out \
    >> CV_phyps2e_m2n1l5.out
done

# Add sample ids to proportions
for k in {1..12}
do
    paste -d " " $base/0_metadata/ids_phyps2e.txt admcv/phyps2e_m2n1l5.${k}.Q |
        sed 's/ $//g' \
        > AdmcvProp_phyps2e_m2n1l5_k${k}.tsv
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


### ============================================================================
### X. Convert VCF to Beagle

#!/bin/bash

#SBATCH --job-name=beagle
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=16G
#SBATCH --time=0-01:00  # D-HH:MM
#SBATCH --output=log/1_beagle_%j.out
#SBATCH --error=log/1_beagle_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-env/8.3
ml VCFtools/0.1.16-GCC-8.3.0

BASE=/gss/work/haex1482/phylo2/6_admix

vcftools \
    --gzvcf $BASE/../2_genotyping/out/9_filt/phyps2_snpsfilt.vcf.gz \
    --chr LG01 \
    --BEAGLE-PL \
    --stdout > $BASE/beagle/phyps2_chr1_beagle.pl

# Manual (= ids_phyps2.txt)
head -n 1 $BASE/beagle/phyps2_chr1_beagle.pl | cut -f 4- | sed 's/\t/\n/g' | uniq > $BASE/samples.txt


### ============================================================================
### Y. Calculate admixture proportions (NGSadmix)

#!/bin/bash

#SBATCH --job-name=ngsadm
#SBATCH --partition=all_nodes.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/2_ngsadm_%j.out
#SBATCH --error=log/2_ngsadm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-env/8.3 HTSlib/1.15.1-GCC-8.3.0

BASE=/gss/work/haex1482/phylo2/6_admix

for k in {2..19}
do
    NGSadmix \
        -likes $BASE/beagle/phyps2_chr1_beagle.pl \
        -K ${k} \
        -o $BASE/ngsadm/phyps2_chr1_k${k} \
        -P 10

    paste -d " " $BASE/samples.txt $BASE/ngsadm/phyps2_chr1_k${k}.qopt |
    sed 's/ $//g' \
    > $BASE/ngsadm/phyps2_chr1_k${k}_ngsadmProp.tsv
done

# Plotted with R/ngsadm_phyps2.R
