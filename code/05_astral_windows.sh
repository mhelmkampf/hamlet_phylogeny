### ============================================================================
### phylo2
### Phylogenetic reconstruction based on genomic-window approach (ASTRAL)
### ============================================================================

### Preparations

mkdir $WORK/phylo2/3_phylogeny/windows/phylo2k

BASE=$WORK/phylo2/3_phylogeny/windows/phylo2k

cd $BASE

mkdir 1_vcf
mkdir 2_fas
mkdir 3_aln
mkdir 4_trees
mkdir 5_astral


### Generate BED file with randomly drawn, non-overlapping windows (REL_COV >= 0.80 & COV_HYP >= 50)

## See select_windows_bed_v2.R (seed 27)
#> win_2000x5kb_s27.bed
#> 10 Mb sequence (1.6% of assembly excluding missing sites)


### ============================================================================
### 1. Extract sites from VCF (loop & array)

#!/bin/bash

#SBATCH --job-name=extract
#SBATCH --partition=carl.p
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/a_extract_%A_%a.out
#SBATCH --error=log/a_extract_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml VCFtools

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

PER_TASK=20

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM ; run<=END_NUM ; run++ ))
do
  echo This is task $SLURM_ARRAY_TASK_ID, run number $run
  chr=$(awk -v line="$run" 'BEGIN { FS = "\t" } ; NR==line+1 { print $1 }' $BASE/win_2000x5kb_s27.bed)
  sta=$(awk -v line="$run" 'BEGIN { FS = "\t" } ; NR==line+1 { print $2 }' $BASE/win_2000x5kb_s27.bed)
  end=$(awk -v line="$run" 'BEGIN { FS = "\t" } ; NR==line+1 { print $3 }' $BASE/win_2000x5kb_s27.bed)
  printf -v i "%04d" $run
  vcftools \
  --gzvcf /gss/work/haex1482/phylo2/2_genotyping/out/8_geno/phylo2_all-indel.vcf.gz \
  --chr "$chr" \
  --from-bp "$sta" \
  --to-bp "$end" \
  --recode --stdout |
  grep -v "##" |
  bgzip > $BASE/1_vcf/window_"$i"_s27.vcf.gz
done

#> use array=1-200 and PER_TASK=10 next time
#> 1960 windows completed (2 arrays failed, see fix below)


### ============================================================================
### 2. Convert to Fasta (loop & array)

#!/bin/bash

#SBATCH --job-name=convert
#SBATCH --partition=carl.p
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/b_convert_%A_%a.out
#SBATCH --error=log/b_convert_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml VCFtools

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

PER_TASK=100

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM ; run<=END_NUM ; run++ ))
do
  echo This is task $SLURM_ARRAY_TASK_ID, run number $run
  printf -v i "%04d" $run
  bgzip -cd $BASE/1_vcf/window_"$i"_s27.vcf.gz |
  vcf-to-tab > $BASE/2_fas/window_"$i"_s27.tab
  perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i $BASE/2_fas/window_"$i"_s27.tab > $BASE/2_fas/window_"$i"_s27.fas
  rm $BASE/2_fas/window_"$i"_s27.tab*
done

# ------------------------------------------------------------------------------
# Fix windows 141-180

#!/bin/bash

#SBATCH --job-name=fix
#SBATCH --partition=carl.p
#SBATCH --array=1-40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/ab_fix_%A_%a.out
#SBATCH --error=log/ab_fix_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml VCFtools

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

PER_TASK=1

START_NUM=$(( ($SLURM_ARRAY_TASK_ID + 139) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK + 140))

echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM ; run<=END_NUM ; run++ ))
do
  echo This is task $SLURM_ARRAY_TASK_ID, run number $run
  chr=$(awk -v line="$run" 'BEGIN { FS = "\t" } ; NR==line+1 { print $1 }' $BASE/win_2000x5kb_s27.bed)
  sta=$(awk -v line="$run" 'BEGIN { FS = "\t" } ; NR==line+1 { print $2 }' $BASE/win_2000x5kb_s27.bed)
  end=$(awk -v line="$run" 'BEGIN { FS = "\t" } ; NR==line+1 { print $3 }' $BASE/win_2000x5kb_s27.bed)
  printf -v i "%04d" $run
  vcftools \
  --gzvcf /gss/work/haex1482/phylo2/2_genotyping/out/8_geno/phylo2_all-indel.vcf.gz \
  --chr "$chr" \
  --from-bp "$sta" \
  --to-bp "$end" \
  --recode --stdout |
  grep -v "##" |
  bgzip > $BASE/1_vcf/window_"$i"_s27.vcf.gz
  bgzip -cd $BASE/1_vcf/window_"$i"_s27.vcf.gz |
  vcf-to-tab > $BASE/2_fas/window_"$i"_s27.tab
  perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i $BASE/2_fas/window_"$i"_s27.tab > $BASE/2_fas/window_"$i"_s27.fas
  rm $BASE/2_fas/window_"$i"_s27.tab*
done


### ============================================================================
### 3. Align sequences (loop & array)

#!/bin/bash

#SBATCH --job-name=c_align
#SBATCH --partition=carl.p
#SBATCH --array=1-40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00  # D-HH:MM
#SBATCH --output=log/c_align_%A_%a.out
#SBATCH --error=log/c_align_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3 MAFFT/7.475-GCC-8.3.0-with-extensions

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

PER_TASK=50

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM ; run<=END_NUM ; run++ ))
do
  echo This is task $SLURM_ARRAY_TASK_ID, run number $run
  printf -v i "%04d" $run
  mafft --auto $BASE/2_fas/window_"$i"_s27.fas > $BASE/3_aln/window_"$i"_s27.aln
done


### ============================================================================
### 4a. Phylogenetic analysis (loop & array)

#!/bin/bash

#SBATCH --job-name=d_iqtree
#SBATCH --partition=carl.p
#SBATCH --array=1-200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00  # D-HH:MM
#SBATCH --output=log/d_iqtree_%A_%a.out
#SBATCH --error=log/d_iqtree_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

PER_TASK=10

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM ; run<=END_NUM ; run++ ))
do
  echo This is task $SLURM_ARRAY_TASK_ID, run number $run
  printf -v i "%04d" $run
  iqtree2 \
    -s $BASE/3_aln/window_"$i"_s27.aln \
    -m TEST \
    --prefix $BASE/4_trees/no/window_"$i"_s27.iqtree \
    -T 1
done


### ============================================================================
### 5. ASTRAL summary tree (loop & array)

#!/bin/bash

#SBATCH --job-name=e_astral
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=5-00:00  # D-HH:MM
#SBATCH --output=log/e_astral_%j.out
#SBATCH --error=log/e_astral_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

cat $BASE/4_trees/no/window_*_s27.iqtree.treefile > $BASE/5_astral/phylo2k_s27_genetrees.tre

java -D"java.library.path=/user/haex1482/apps/ASTRAL-MP_v5.15.5/lib" \
  -jar ~/apps/ASTRAL-MP_v5.15.5/astral.5.15.5.jar \
  -i $BASE/5_astral/phylo2k_s27_genetrees.tre \
  -o $BASE/5_astral/phylo2k_s27_astral.tre \
  2> $BASE/5_astral/phylo2k_s27_astral.log \
  -T 24


### ============================================================================
### 4b. Phylogenetic analysis with bootstrapping (loop & array)

#!/bin/bash

#SBATCH --job-name=f_iqtree
#SBATCH --partition=carl.p
#SBATCH --array=1-200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5-00:00  # D-HH:MM
#SBATCH --output=log/f_iqtree_%A_%a.out
#SBATCH --error=log/f_iqtree_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

PER_TASK=10

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

for (( run=$START_NUM ; run<=END_NUM ; run++ ))
do
  echo This is task $SLURM_ARRAY_TASK_ID, run number $run
  printf -v i "%04d" $run
  iqtree2 \
    -s $BASE/3_aln/window_"$i"_s27.aln \
    -m TEST \
    -B 1000 \
    --prefix $BASE/4_trees/ufb/window_"$i"_s27_ufb.iqtree \
    -T 1
done


# ------------------------------------------------------------------------------
# Calculate mean support value per window
for i in *.contree
do
  echo $i $( \
  grep -oE "[0-9]+:" $i |
  sed 's/://' |
  awk '{ sum += $1 ; next } END { print sum / NR }') \
  >> phylo2k_s27_meanufb.tsv
done

# Concatenate gene trees meeting average posterior probability threshold
cat $(awk '$2 >= 30 { print $1 }' phylo2k_s27_meanufb.tsv |
      sed 's/contree/treefile/') \
    >> ../../5_astral/phylo2k_s27_genetrees_minufb30.tre
#> 1805 gene trees

cat $(awk '$2 >= 50 { print $1 }' phylo2k_s27_meanufb.tsv |
      sed 's/contree/treefile/') \
    >> ../../5_astral/phylo2k_s27_genetrees_minufb50.tre
#> 1321 gene trees

# Top 10 gene trees by PP support
sort -nr -k2 phylo2k_s27_meanufb.tsv | head -n 3
# window_0439_s27_pp.iqtree.contree 87.3494   # 1224 parsimony-informative sites
# window_0381_s27_pp.iqtree.contree 86.7319   # 1121 parsimony-informative sites
# window_0277_s27_pp.iqtree.contree 86.6295   # 1125 parsimony-informative sites

sort -n -k2 phylo2k_s27_meanufb.tsv | head -n 3
# window_1619_s27_pp.iqtree.contree 15.6655   # 519 parsimony-informative sites
# window_0339_s27_pp.iqtree.contree 15.9462   # 528 parsimony-informative sites
# window_0813_s27_pp.iqtree.contree 17.2837   # 620 parsimony-informative sites


### ============================================================================
### 5b. ASTRAL summary tree based on bootstrap threshold (loop & array)

#!/bin/bash

#SBATCH --job-name=g_astral
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=3-00:00  # D-HH:MM
#SBATCH --output=log/g_astral_%j.out
#SBATCH --error=log/g_astral_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

BASE=/gss/work/haex1482/phylo2/3_phylogeny/windows/phylo2k

java -D"java.library.path=/user/haex1482/apps/ASTRAL-MP_v5.15.5/lib" \
  -jar ~/apps/ASTRAL-MP_v5.15.5/astral.5.15.5.jar \
  -i $BASE/5_astral/phylo2k_s27_genetrees_minufb30.tre \
  -o $BASE/5_astral/phylo2k_s27_astral_minufb30.tre \
  2> $BASE/5_astral/phylo2k_s27_astral_minufb30.log \
  -T 24

java -D"java.library.path=/user/haex1482/apps/ASTRAL-MP_v5.15.5/lib" \
  -jar ~/apps/ASTRAL-MP_v5.15.5/astral.5.15.5.jar \
  -i $BASE/5_astral/phylo2k_s27_genetrees_minufb50.tre \
  -o $BASE/5_astral/phylo2k_s27_astral_minufb50.tre \
  2> $BASE/5_astral/phylo2k_s27_astral_minufb50.log \
  -T 24


# ------------------------------------------------------------------------------
# Mean nodal support (quadripartition PP) in ASTRAL trees:
for i in *astral*.tre
do
  echo $i $( \
  grep -oE "[0-9]+:" $i |
  sed 's/://' |
  awk '{ sum += $1 ; next } END { print sum / NR }')
done

#> phylo2k_s27_astral.tre            44.9729
#> phylo2k_s27_astral_minufb30.tre   43.5482
#> phylo2k_s27_astral_minufb50.tre   43.7560

# Mean branch length (coalescent units) in ASTRAL trees:
for i in *astral*.tre
do
  echo $i $( \
  grep -oE ":[0-9/.]+" $i |
  sed 's/://' |
  awk '{ sum += $1 ; next } END { print sum / NR }')
done

#> phylo2k_s27_astral.tre            0.344106
#> phylo2k_s27_astral_minufb30.tre   0.379738
#> phylo2k_s27_astral_minufb50.tre   0.310098

## Repeat with branch length (concordance)?

# Note: ufb originally named pp