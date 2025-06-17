### ================================================================================
### Radiation with reproductive isolation in the near-absence of phylogenetic signal
### 04b. Species tree inference (SNPs, SVDQuartets)
### By Martin Helmkampf, last edited 2025-06-16
### ================================================================================

### Preparations

base=$WORK/phylo2

mkdir $base/3_phylogeny/svdq

cd $base/3_phylogeny/svdq


### ============================================================================
### 1. Convert VCF to Nexus file

gzip -cd $base/2_genotyping/out/9_filt/phylo-snp_5kb.vcf.gz \
    > phylo-snp_5kb.vcf

ruby ~/apps/convert_vcf_to_nexus.rb phylo-snp_5kb.vcf phylo-snp_5kb.nex

# rm phylo-snp_5kb.vcf

# Obtain outgroup taxon numbers
grep '#CHROM' phylo-snp_5kb.vcf |
    sed 's/\s/\n/g' |
    tail -n +10 \
    > sample_numbers.txt

grep 'tor\|tab\|tig' sample_numbers.txt
# 135	20478tabhon
# 136	20480tabhon
# 137	20481tighon
# 210	28393torboc
# 249	Bocas16.3torboc
# 250	Bocas16.4torboc
# 290	PL17_21tigboc
# 335	s_tort_3torboc


### ============================================================================
### 2. Run PAUP

#!/bin/bash

#SBATCH --job-name=svdq1
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=03-00   # DD-HH
#SBATCH --output=sl_%j_svdq1.out
#SBATCH --error=sl_%j_svdq1.err

paup -n <<EOF
log file=phylo-snp_5kb_svdq_1K.log;
execute phylo-snp_5kb.nex;
outgroup 135 136 137 210 249 250 290 335;
svdq nquartets=490000000 ambigs=distribute seed=444 nthreads=128;
savetrees file=phylo-snp_5kb_svdq_1K.nex;
quit;
EOF


### ----------------------------------------------------------------------------
### Bootstrap

#!/bin/bash

#SBATCH --job-name=svdqbs
#SBATCH --partition=rosa.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=0
#SBATCH --time=05-00   # DD-HH
#SBATCH --output=sl_%j_svdqbs.out
#SBATCH --error=sl_%j_svdqbs.err

paup -n <<EOF
log file=phylo-snp_5kb_svdq_bs100C.log;
execute phylo-snp_5kb.nex;
outgroup 135 136 137 210 249 250 290 335;
svdq nquartets=5000000 ambigs=distribute bootstrap nreps=200 seed=999 nthreads=128 treefile=phylo-snp_5kb_svdq_bs200C.nex;
quit;
EOF

# 100B: nquartets=100000 seed=456 nreps=100
# 200C: nquartets=5000000 bootstrap nreps=200 seed=999 > better resolution


### ----------------------------------------------------------------------------
### Consensus tree (executed on Paup command line)

getTrees file=phylo-snp_5kb_svdq_bs200C.nex;   # maxtrees=200
contree /majrule file=phylo-snp_5kb_svdq_bs200C_maj.nex;

# Draw support values on tree: https://ib.berkeley.edu/courses/ib200/2018/labs/04/lab04.pdf ?


### ============================================================================
### 3. quartetsampling (QS)

mkdir qs && cd qs

python3 /user/haex1482/apps/vcf2phylip/vcf2phylip.py -i ../phylo-snp_5kb.vcf
mv ../phylo-snp_5kb.min4.phy phylo-snp_5kb.phy

# tree phylo-snp_5kb_svdq_1K.nex converted to Newick format by write.tree()
sed 's/NA//g' phylo-snp_5kb_svdq_1K.nwk > phylo-snp_5kb_svdq_1K.ed.nwk

conda activate python3.7

ml RAxML-NG/1.2.0-GCC-13.1.0

python3 /user/haex1482/apps/quartetsampling/pysrc/quartet_sampling.py \
    --tree phylo-snp_5kb_svdq_1K.ed.nwk \
    --align phylo-snp_5kb.phy \
    --reps 100 \
    --lnlike 2 \
    --threads 24

rename RESULT phylo-snp_5kb_svdq_1K_qs *
