### ============================================================================
### phylo2
### Phasing (Shapeit)
### ============================================================================

### Process overview

# - make bamlist file (bash)
# - split by LG (vcftools)
# - extract PIRs (Shapeit)
# - phase (Shapeit)
# - merge LGs (vcftools)


### ============================================================================
### Preparations

base=$WORK/phylo2

mkdir $base/4_phasing

cd $base/4_phasing

mkdir 1_split
mkdir 2_pirs
mkdir 3_phase
mkdir 4_merge
mkdir log
mkdir scripts

# unaligned BAM files in --
# deduplicated BAM files in $DATA/projects/5_Phylo/2_genotyping/out/4_dedup

## -----------------------------------------------------------------------------
## Create bamlist

cd $DATA/projects/5_Phylo/2_genotyping/out/4_dedup

for f in *.bam
do
    DEDUP=$DATA/projects/5_Phylo/2_genotyping/out/4_dedup/${f}
    SAMPLE=${f%_dedup.bam}
    # UBAM=$WORK/phylo2/2_genotyping/out/1_ubam/${SAMPLE}_ubam.bam
    printf $SAMPLE" "$DEDUP"\n"
done > $base/bamlist_phylo2.tsv
cd -

# *** Problem: dedup files are merged (342), ubam files are not (552) ***

for i in {01..24}
do
    awk -v lg=$i '{ print $0, "LG"lg }' $base/bamlist_phylo2.tsv
done |
sort -k 3 > $base/bamlist_LGs_phylo2.tsv


### ============================================================================
### 1. Split into chromosomes (VCFtools)

#!/bin/bash

#SBATCH --job-name=1_split
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:00  # D-HH:MM
#SBATCH --output=log/1_split_%A_%a.out
#SBATCH --error=log/1_split_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-env/8.3 VCFtools/0.1.16-GCC-8.3.0

base=/gss/work/haex1482/phylo2/4_phasing
SEQ=($(seq -w 1 24))

vcftools \
    --gzvcf $WORK/phylo2/2_genotyping/out/9_filt/phylo2_snpsfilt.vcf.gz \
    --chr LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]} \
    --recode \
    --stdout |
bgzip > $base /1_split/phylo2_snpsfilt_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz
tabix -p vcf $base /1_split/phylo2_snpsfilt_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz


### ============================================================================
### 2. Extract PIRs (extractPIRs)

#!/bin/bash

#SBATCH --job-name=2_extract
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-02:00  # D-HH:MM
#SBATCH --output=log/2_extract_%A_%a.out
#SBATCH --error=log/2_extract_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

base=/gss/work/haex1482/phylo2/4_phasing
SEQ=($(seq -w 1 24))

extractPIRs \
    --bam $base/bamlist_LGs_phylo2.tsv \
    --vcf $base /1_split/phylo2_snpsfilt_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    --out $base /2_pirs/pirlist_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.txt \
    --base-quality 20 \
    --read-quality 15 \
    --threads 12


### ============================================================================
### 3. Phase (Shapeit v2)

#!/bin/bash

#SBATCH --job-name=3_phase
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-00:00  # D-HH:MM
#SBATCH --output=log/3_phase_%A_%a.out
#SBATCH --error=log/3_phase_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

base=/gss/work/haex1482/phylo2/4_phasing
SEQ=($(seq -w 1 24))

shapeit -assemble \
    --input-vcf $base /1_split/phylo2_snpsfilt_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    --input-pir $base /2_pirs/pirlist_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.txt \
    --output-max $base /3_phase/phased_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]} \
    --force \
    --thread 12
   #--output-log $base /3_phase/phased_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}_asm.log \

shapeit -convert \
    --input-hap $base /3_phase/phased_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]} \
    --output-vcf $base /3_phase/phased_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    --thread 12
 #--output-log $base /3_phase/phased_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}_con.log \

# rm $base/*.log


### ============================================================================
### 4. Merge (VCFtools)

#!/bin/bash

#SBATCH --job-name=4_merge
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-06:00  # D-HH:MM
#SBATCH --output=log/4_merge_%j.out
#SBATCH --error=log/4_merge_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-env/8.3 VCFtools/0.1.16-GCC-8.3.0

base=/gss/work/haex1482/phylo2/4_phasing

vcf-concat $base /3_phase/phased_LG*.vcf.gz |
grep -v ^\$ |
tee $base /4_merge/phylo2_phased.vcf |
vcftools --vcf - --mac 2 --recode --stdout > $base /4_merge/phylo2_phased_mac2.vcf

bgzip $base /4_merge/phylo2_phased.vcf.gz
bgzip $base /4_merge/phylo2_phased_mac2.vcf.gz

tabix -p vcf $base /4_merge/phylo2_phased.vcf.gz
tabix -p vcf $base /4_merge/phylo2_phased_mac2.vcf.gz


## -----------------------------------------------------------------------------
## Relabel samples

cd $base

bcftools reheader \
    -s relabel.txt \
    --threads 8 \
    -o phylo2e_phased.vcf.gz \
    phylo2_phased.vcf.gz

tabix -p vcf phylo2e_phased.vcf.gz



### ============================================================================
### Whatshap / Shapeit4 approach (uses recombination map) *** not executed ***
### ============================================================================

### 2b. Whatshap

#!/bin/bash

#SBATCH --job-name=2b_whap
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250G
#SBATCH --time=5-00:00  # D-HH:MM
#SBATCH --output=log/2b_whap_%A_%a.out
#SBATCH --error=log/2b_whap_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

base=/gss/work/haex1482/phylo2/4_phasing
SEQ=($(seq -w 1 24))

whatshap phase \
  $base /1_split/phylo2_snpsfilt_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
  $DATA/projects/5_Phylo/2_genotyping/out/4_dedup/*.bam \
  --reference=$base/../2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
  -o $base/whap/phylo2_whap_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
  --tag=PS


### ============================================================================
### 3b. Shapeit v4

#!/bin/bash

#SBATCH --job-name=3b_phase
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=5-00:00  # D-HH:MM
#SBATCH --output=log/3b_phase_%A_%a.out
#SBATCH --error=log/3b_phase_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml hpc-env/8.3 SHAPEIT4/4.2.2-foss-2019b

base=/gss/work/haex1482/phylo2/4_phasing

shapeit4.2 \
  --input $base /1_split/phylo2_snpsfilt_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
  --output $base/phase4/phased4_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
  --map $base/Hpue_map1_prelim.gmap \
  --region LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]} \
  --use-PS 0.0001 \
  --sequencing \
  --thread 12
