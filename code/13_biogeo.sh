### ============================================================================
### phylo2
### Biogeographic analysis (BioGeoBEARS)
### ============================================================================

### Preparations

base=$WORK/phylo2

cd $base/3_phylogeny/iqtree


### ============================================================================
### 1. Dating (LSD2 method)

# Dating file
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
