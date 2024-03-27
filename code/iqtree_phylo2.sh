### ============================================================================
### phylo2
### Phylogenetic reconstruction based on concatenated SNPs (IQ-TREE)
### ============================================================================

### Preparations

base=$WORK/phylo2

mkdir $base/3_phylogeny/iqtree

cd $base/3_phylogeny/iqtree


### ============================================================================
### 1. Convert VCF to Fasta

#!/bin/bash

#SBATCH --job-name=fasta
#SBATCH --partition=rosa_express.p
#SBATCH --nodes=1
#SBATCH --time=00-02   # DD-HH
#SBATCH --output=/dev/null
#SBATCH --error=sl_%j_fasta.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml VCFtools/0.1.16-GCC-13.1.0

base=$WORK/phylo2

gzip -cd $base/2_genotyping/out/9_filt/phylo2e_m2k5.vcf.gz \
    > phylo2e_m2k5.vcf

vcf-to-tab < \
    phylo2e_m2k5.vcf |
    sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' \
    > phylo2e_m2k5.tab

perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl \
    -i phylo2e_m2k5.tab \
    > phylo2e_m2k5.fas

# Confirm sequences are of equal length (aligned)
bioawk -c fastx '{ print $name, length($seq) }' phylo2e_m2k5.fas   # 110436 aligned sites (v1: 109796)

rm phylo2e_m2k5.vcf
rm phylo2e_m2k5.tab*


### ============================================================================
### 2. Run IQ-TREE

ml IQ-TREE/2.2.2.7-gompi-2023a

iqtree2 \
    -s phylo2e_m2k5.fas \
    -m GTR+ASC

rm *.ckp.gz

#> +ASC failed due to 27358 invariable sites present (see phylo2e_m2k5.fas.log)
#> Variable sites printed to phylo2e_m2k5.fas.varsites.phy
#> 70505 parsimony-informative sites

#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH --partition=rosa.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00-06   # DD-HH
#SBATCH --output=sl_%j_iqtree.out
#SBATCH --error=sl_%j_iqtree.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml IQ-TREE/2.2.2.7-gompi-2023a

iqtree2 \
    -s phylo2e_m2k5.fas.varsites.phy \
    -m GTR+ASC \
    -allnni \
    -nbest 10 \
    -B 1000 \
    -bnni \
    -T AUTO \
    --prefix iqtree_phylo2e_m2k5_GTRL_1A


### ============================================================================
### 3. Dating (LSD method)

# Dating (v1 with tabs)
cat > dates_ftol.txt <<EOF
28393torboc,20480tabhon -10.68
28393torboc,62559floarc -25.93
28393torboc,PL17_21tigboc -40.96
EOF

#!/bin/bash

#SBATCH --job-name=lsd
#SBATCH --partition=rosa.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03-00   # DD-HH
#SBATCH --output=sl_%j_lsd.out
#SBATCH --error=sl_%j_lsd.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

# ml IQ-TREE/2.2.2.7-gompi-2023a (run with build v2.2.5)

iqtree2 \
    -s phylo2e_m2k5.fas.varsites.phy \
    -te iqtree_phylo2e_m2k5_GTRL_1A.treefile \
    -o "PL17_21tigboc,20481tighon" \
    -m GTR+ASC \
    -nt AUTO \
    --date dates_ftol.txt \
    --date-tip 0 \
    --date-ci 100 \
    --date-options "-u 0" \
    --prefix lsd_phylo2e_m2k5_GTRL_1A
