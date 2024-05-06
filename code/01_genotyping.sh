### ============================================================================
### phylo2
### Data compilation and genotyping (GATK)
### ============================================================================

### Preparations

mkdir $WORK/phylo2/0_metadata
mkdir $WORK/phylo2/1_rawdata
mkdir $WORK/phylo2/2_genotyping

BASE=$WORK/phylo2/2_genotyping

cd $BASE

## Data summary (n = 337 individuals, 342 datasets, 552 lanes, 1104 files)
## wgs17:            232 samples, mostly 2 files each              *** low coverage: L20616, PL17_79, PL17_98, PL17_101 ***
## phylo2:           104 samples, 6 files each (624 files total)   *** S135 dropped ***
## betancur:           6 samples, 2 files each  (12 files total)
##
## Species / location breakdown see R/tally_phylo2.R, producing phylo2_tally_12.pdf / .tsv,
## with new seqdata containing updated species / location tags in 1_Resources/2_Metadata/seqdata_new_wgs18.tsv
##
## Samples with duplicate libraries: 19519affgun, 28393torboc, PL17_79abepri, PL17_98indbel, PL17_111indbel


## Compile raw data (symlinks, in $BASE/1_rawdata)

# wgs17 samples were chosen in R/choose_wgs17_phylo2.R, including:

# all in $DATA/shared/1_rawdata/wgs_12 (107 individuals, 218 files -- 2 in duplicate)
ln -s /nfs/data/haex1482/shared/1_rawdata/wgs_12/* .

# select samples in $DATA/shared/1_rawdata/wgs_3 (48 individuals, 96 files)
while read SAMPLE
do
  ln -s /nfs/data/haex1482/shared/1_rawdata/wgs_3/${SAMPLE}* .
done < wgs_3.ids

# all in $DATA/shared/1_rawdata/wgs_456 (21 individuals, 42 files)
ln -s /nfs/data/haex1482/shared/1_rawdata/wgs_456/* .

# select samples in /nfs/data/doau0129/GxP (56 individuals, 112 files)
while read SAMPLE
do
  ln -s /nfs/data/doau0129/GxP/${SAMPLE}* .
done < wgs_gxp.ids

# phylo2 in /nfs/data/doau0129/phylo2/phylo2_raw (104 individuals, 624 files)
ln -s /nfs/data/doau0129/phylo2/phylo2_raw/*.gz .
rm *S135*

# Additional samples by Betancur downloaded from ENA and SRA (6 individuals, 12 files) *** 7 individuals?
# FL0835 SRR17839752 pro san
# FL0331 SRR18184334 chl pri
# FL0880 SRR17839729 ran san
# FL0836 SRR18184176 gut san
# FL0839 SRR19070321 abe san (raw files renamed)
# FL0318 SRR19070349 ind pri (raw files renamed)
ln -s /user/haex1482/data/shared/1_rawdata/betancur/* .


## Compile metadata
find ../1_rawdata/*.f*.gz > fofn/0_raw.fofn
# metadata file for genotyping prepared with R/prepare_meta_phylo2.R,
# based on seqdata_new_wgs18.tsv, betancur_ids.csv, and 0_raw.fofn


## Unzip and index reference genome assembly
ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8
ml SAMtools/1.9-GCC-8.3.0
ml BWA/0.7.17-GCC-8.3.0

gzip -dc /nfs/data/haex1482/shared/2_refdata/Hpue_genome_unmasked_01.fas.gz > $BASE/ref/Hpue_genome_unmasked_01.fasta
bwa index $BASE/ref/Hpue_genome_unmasked_01.fasta
gatk CreateSequenceDictionary -R $BASE/ref/Hpue_genome_unmasked_01.fasta
samtools faidx $BASE/ref/Hpue_genome_unmasked_01.fasta


## Set up directories for genotyping
for i in log out
do
  mkdir $i
  for j in 1_ubam 2_adap 3_map 4_dedup 5_cov 6_like 7_coho 8_geno 9_filt
  do
    mkdir $i/$j
  done
done


### ============================================================================
### 1. Convert Fastq to unaligned BAM

#!/bin/bash

#SBATCH --job-name=1_ubam
#SBATCH --partition=carl.p
#SBATCH --array=1-552
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-01:30  # D-HH:MM
#SBATCH --output=log/1_ubam/1_ubam_%A_%a.out
#SBATCH --error=log/1_ubam/1_ubam_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
META=$BASE/../0_metadata/metadata_phylo2.csv
LINES=$(head $META -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS="," read FWD REV SAMPLE READGROUP LIBRARY PLATFORM MODEL <<< $LINES

gatk --java-options "-Xmx20G" \
    FastqToSam \
    -F1 $BASE/../1_rawdata/$FWD \
    -F2 $BASE/../1_rawdata/$REV \
    -O $BASE/out/1_ubam/${READGROUP}_ubam.bam \
    --SAMPLE_NAME $SAMPLE \
    --READ_GROUP_NAME $READGROUP \
    --LIBRARY_NAME $LIBRARY \
    --PLATFORM $PLATFORM \
    --PLATFORM_MODEL $MODEL \
    --TMP_DIR $BASE/tmp

## Manually create list of files
cd $BASE/out/1_ubam ; find *_ubam.bam > $BASE/fofn/1_ubam.fofn


### ============================================================================
### 2. Mark adapters

#!/bin/bash

#SBATCH --job-name=2_adap
#SBATCH --partition=carl.p
#SBATCH --array=1-552
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-01:30  # D-HH:MM
#SBATCH --output=log/2_adap/2_adap_%A_%a.out
#SBATCH --error=log/2_adap/2_adap_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
FOFN=$BASE/fofn/1_ubam.fofn
UBAM=$(head $FOFN -n $SLURM_ARRAY_TASK_ID | tail -n 1)
READGROUP=${UBAM%_ubam.bam}

gatk --java-options "-Xmx20G" \
   MarkIlluminaAdapters \
   -I $BASE/out/1_ubam/$UBAM \
   -O $BASE/out/2_adap/${READGROUP}_adap.bam \
   -M $BASE/out/2_adap/${READGROUP}_adap.txt \
   --TMP_DIR $BASE/tmp

## Manually create list of files
cd $BASE/out/2_adap ; find *_adap.bam > $BASE/fofn/2_adap.fofn


### ============================================================================
### 3. Map to reference genome and sort

#!/bin/bash

#SBATCH --job-name=3_map
#SBATCH --partition=carl.p
#SBATCH --array=1-552
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=3-00:00  # D-HH:MM
#SBATCH --output=log/3_map/3_map_%A_%a.out
#SBATCH --error=log/3_map/3_map_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8
ml BWA/0.7.17-GCC-8.3.0

BASE=/gss/work/haex1482/phylo2/2_genotyping
FOFN=$BASE/fofn/2_adap.fofn
ADAP=$(head $FOFN -n $SLURM_ARRAY_TASK_ID | tail -n 1)
READGROUP=${ADAP%_adap.bam}

gatk --java-options "-Xmx60G" \
    SamToFastq \
    -I $BASE/out/2_adap/$ADAP \
    -F /dev/stdout \
    --INTERLEAVE true \
    --INCLUDE_NON_PF_READS true \
    --TMP_DIR $BASE/tmp |

bwa mem -M -t 8 \
    -p $BASE/ref/Hpue_genome_unmasked_01.fasta /dev/stdin |

gatk --java-options "-Xmx60G" \
    MergeBamAlignment \
    -ALIGNED /dev/stdin \
    -UNMAPPED $BASE/out/1_ubam/${READGROUP}_ubam.bam \
    -O $BASE/out/3_map/${READGROUP}_map.bam \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --ADD_MATE_CIGAR true \
    --ALIGNED_READS_ONLY false \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --ATTRIBUTES_TO_RETAIN X0 \
    --CLIP_ADAPTERS false \
    --EXPECTED_ORIENTATIONS FR \
    --MAX_RECORDS_IN_RAM 2000000 \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --SORT_ORDER "coordinate" \
    --UNMAP_CONTAMINANT_READS true \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --VALIDATION_STRINGENCY SILENT \
    --TMP_DIR $BASE/tmp

## Manually create list of files
cd $BASE/out/3_map ; find *_map.bam > $BASE/fofn/3_map.fofn


## -----------------------------------------------------------------------------
## Prepare file with lanes to merge (involves manual steps)

awk 'BEGIN { FS = "," } ;
    { if (!seen[$5]++) print $5","$4 ; else print $4 }' metadata_phylo2.csv | sort > lane_tmp.csv

# File prepared manually *** automate in R ***
# $ >> _map.bam (remove from last line)
# ([\w-]+),([\w.]+)\n([\w.]+)\n([\w.]+) >> $1,$2_map.bam $3_map.bam $4_map.bam
# fix 2 outliers
# add -I
#> lane_merge.csv

# Regex solution in sed / perl (in progress)
sed -E 's/(.+), (.+)\n(.+)\n(.+)/\1, \2 \3 \4/g' tmp.csv

sed -E 's/([0-9a-z-]+),\.+/\1/g' tmp.csv

perl -pe 's/([\w-]+), ([\w.]+)\n([\w.]+)\n([\w.]+)/$1, $2 $3 $4/g' tmp.csv

# Complete solution in awk (in progress)
awk 'BEGIN { FS = "," } ;
    { (!seen[$5]) ? a = $5", "$4 : a = $4 ; print a }' metadata_phylo2.csv

# Solution in R

#> Lanes merged to 342 libraries


### ============================================================================
### 4. Sort and mark duplicates, merge lanes (*** may fail for some samples ***)

#!/bin/bash

#SBATCH --job-name=4_dedup
#SBATCH --partition=carl.p
#SBATCH --array=1-342
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/4_dedup/4_dedup_%A_%a.out
#SBATCH --error=log/4_dedup/4_dedup_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
LIST=$BASE/../0_metadata/lane_merge.csv
LINES=$(head $LIST -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS="," read LIBRARY LANES <<< $LINES

cd $BASE/out/3_map/

gatk --java-options "-Xmx40G" \
    MarkDuplicates \
    $LANES \
    -O $BASE/out/4_dedup/${LIBRARY}_dedup.bam \
    -M $BASE/out/4_dedup/${LIBRARY}_dedup.txt \
    --REMOVE_DUPLICATES false \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
    --TMP_DIR $BASE/tmp

gatk --java-options "-Xmx40G" \
    BuildBamIndex \
    -I $BASE/out/4_dedup/${LIBRARY}_dedup.bam

## Manually create list of files
cd $BASE/out/4_dedup ; find *_dedup.bam > $BASE/fofn/4_dedup.fofn

#> 20429ranhon-1 and 20558puehon-1 initially truncated, both rerun separately
#> 18152puebel-1 failed initially, but passed on second run


### ============================================================================
### 5. Calculate coverage

#!/bin/bash

#SBATCH --job-name=5_cov
#SBATCH --partition=carl.p
#SBATCH --array=1-342
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-01:30  # D-HH:MM
#SBATCH --output=log/5_cov/5_cov_%A_%a.out
#SBATCH --error=log/5_cov/5_cov_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
FOFN=$BASE/fofn/4_dedup.fofn
DEDUP=$(head $FOFN -n $SLURM_ARRAY_TASK_ID | tail -n 1)
LIBRARY=${DEDUP%_dedup.bam}

gatk --java-options "-Xmx20G" \
  CollectWgsMetrics \
  -I $BASE/out/4_dedup/$DEDUP \
  -O $BASE/out/5_cov/${LIBRARY}_cov.txt \
  -R $BASE/ref/Hpue_genome_unmasked_01.fasta


## -----------------------------------------------------------------------------
## Prepare coverage table (execute manually)
cd $BASE/out/5_cov
for f in *.txt
do
  lib=${f%_cov.txt}
  cov=$(awk 'FNR == 8 { print $2 }' $f)
  echo $lib $cov >> coverage_phylo2.csv
done

# Plot coverage in R (coverage_phylo2.R)

# Identify libraries not meeting 10x cutoff
awk '($2 < 10) { print $1 }' coverage_phylo2.csv
#> 18436nigboc-1
#> 20616gumhon-1
#> PL17_79abepri-1
#> PL17_98indbel-1

# Pass on libraries for further analysis
awk '($2 >= 10) { print $1"_dedup.bam" }' coverage_phylo2.csv > $BASE/fofn/5_cov.fofn
#> 338 libraries (remaining duplicates: 19519affgun, 28393torboc, PL17_111indbel)

# Mean coverage (before filtering):
awk '{ sum += $2 } END { print sum / NR }' coverage_phylo2.csv
#> 22.1 (after filtering: 22.3)


### ============================================================================
### 6. Calculate haplotype likelihoods

#!/bin/bash

#SBATCH --job-name=6_like
#SBATCH --partition=carl.p
#SBATCH --array=1-338
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/6_like/6_like_%A_%a.out
#SBATCH --error=log/6_like/6_like_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
FOFN=$BASE/fofn/5_cov.fofn
DEDUP=$(head $FOFN -n $SLURM_ARRAY_TASK_ID | tail -n 1)
LIBRARY=${DEDUP%_dedup.bam}

gatk --java-options "-Xmx40G" \
    HaplotypeCaller \
    -I $BASE/out/4_dedup/$DEDUP \
    -O $BASE/out/6_like/${LIBRARY}_g.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    -ERC GVCF

## Manually create list of files
cd $BASE/out/6_like ; find *_g.vcf.gz > $BASE/fofn/6_like.fofn


### ============================================================================
### 7. Combine individual GVCF files into cohort GVCF -- phylo2

#!/bin/bash

#SBATCH --job-name=7_coho
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/7_coho/7_coho_%j.out
#SBATCH --error=log/7_coho/7_coho_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
FOFN=$BASE/fofn/6_like.fofn
INPUT=$(cat $FOFN | sed -e 's/^/-V /g')

cd $BASE/out/6_like

gatk --java-options "-Xmx120G" \
    CombineGVCFs \
    $INPUT \
    -O $BASE/out/7_coho/phylo2_cohort_g.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta

## Clean up cohort VCF *** untested modification ***
# cd $BASE/out/7_coho
# cp phylo2_cohort_g.vcf.gz phylo2_cohort_tmp.vcf.gz
#
# zcat < phylo2_cohort_tmp.vcf.gz |
# grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
# bgzip > phylo2_cohort_g.vcf.gz
#
# rm phylo2_cohort_tmp.vcf.gz


### ============================================================================
### 8.1. Genotype jointly (all sites and SNPs only) -- phylo2

#!/bin/bash

#SBATCH --job-name=8.1_geno
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/8_geno/8.1_geno_%A_%a.out
#SBATCH --error=log/8_geno/8.1_geno_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
SEQ=($(seq -w 1 24))

gatk --java-options "-Xmx120G" \
    GenotypeGVCFs \
    -V $BASE/out/7_coho/phylo2_cohort_g.vcf.gz \
    -L LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]} \
    -O $BASE/out/8_geno/phylo2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --include-non-variant-sites true

gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/8_geno/phylo2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -O $BASE/out/8_geno/phylo2_all-indel_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --select-type-to-exclude INDEL

gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/8_geno/phylo2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -O $BASE/out/8_geno/phylo2_rawsnps_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --select-type-to-include SNP


### ============================================================================
### 8.2. Gather LGs and extract variants stats -- phylo2

#!/bin/bash

#SBATCH --job-name=8.2_geno
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=00-02:00  # D-HH:MM
#SBATCH --output=log/8_geno/8.2_geno_%j.out
#SBATCH --error=log/8_geno/8.2_geno_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
RAW=$(find $BASE/out/8_geno/phylo2_rawsnps_LG*.vcf.gz | awk '{ print "-I", $1 }')
ALL=$(find $BASE/out/8_geno/phylo2_all-indel_LG*.vcf.gz | awk '{ print "-I", $1 }')

# All sites (except indels)
gatk --java-options "-Xmx120g" \
    GatherVcfs \
    $ALL \
    -O $BASE/out/8_geno/phylo2_all-indel.vcf.gz

tabix -p vcf $BASE/out/8_geno/phylo2_all-indel.vcf.gz

# SNPs only (had to be executed manually)
gatk --java-options "-Xmx120g" \
    GatherVcfs \
    $RAW \
    -O $BASE/out/8_geno/phylo2_rawsnps.vcf.gz

tabix -p vcf $BASE/out/8_geno/phylo2_rawsnps.vcf.gz

gatk --java-options "-Xmx120G" \
    VariantsToTable \
    -V $BASE/out/8_geno/phylo2_rawsnps.vcf.gz \
    -O $BASE/stats/genoqual_phylo2_rawsnps.tsv \
    -F CHROM -F POS -F MQ -F QD -F FS -F MQRankSum -F ReadPosRankSum \
    --show-filtered


## -----------------------------------------------------------------------------
## Remove intermediary files (execute manually)
rm $BASE/out/8_geno/phylo2_intermed_*.vcf*
rm $BASE/out/8_geno/phylo2_*LG*


## -----------------------------------------------------------------------------
## Plot distribution of variant quality stats in R (execute manually)

# extract 5% of lines
shuf -n 3500000 genoqual_phylo2_rawsnps.tsv > genoqual_phylo2_rawsnps_3.5M.tsv

# plotting script: qual_distribution_phylo2.R


### ============================================================================
### 9. Filter variants -- phylo2

#!/bin/bash

#SBATCH --job-name=9_filter
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/9_filt/9_filt_%j.out
#SBATCH --error=log/9_filt/9_filt_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8
ml VCFtools/0.1.16-GCC-8.3.0

BASE=/gss/work/haex1482/phylo2/2_genotyping
NOSAM=$(find $BASE/out/6_like/*.vcf.gz | wc -l)
MISS=$((NOSAM / 10))

gatk --java-options "-Xmx120G" \
    VariantFiltration \
    -V $BASE/out/8_geno/phylo2_rawsnps.vcf.gz \
    -O $BASE/out/9_filt/phylo2_snpsint.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --filter-expression "MQ < 57.5" \
    --filter-name "filter_MQ" \
    --filter-expression "QD < 3.0" \
    --filter-name "filter_QD" \
    --filter-expression "FS > 60.0" \
    --filter-name "filter_FS" \
    --filter-expression "MQRankSum < -0.2 || MQRankSum > 0.2" \
    --filter-name "filter_MQRankSum" \
    --filter-expression "ReadPosRankSum < -2.25 || ReadPosRankSum > 2.25 " \
    --filter-name "filter_ReadPosRankSum"

gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/9_filt/phylo2_snpsint.vcf.gz \
    -O $BASE/out/9_filt/phylo2_snpsint2.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --exclude-filtered

vcftools \
    --gzvcf $BASE/out/9_filt/phylo2_snpsint2.vcf.gz \
    --max-missing-count $MISS \
    --max-alleles 2 \
    --stdout \
    --recode |
bgzip > $BASE/out/9_filt/phylo2_snpsfilt.vcf.gz

tabix -p vcf $BASE/out/9_filt/phylo2_snpsfilt.vcf.gz


## -----------------------------------------------------------------------------
## Remove intermediary files (execute manually)
rm $BASE/out/9_filt/*snpsint*


### ============================================================================
### Cleanup and stats

rm -rf tmp
rm 1_ubam/*.bam
rm 2_adap/*.bam

mkdir scripts && mv *.js scripts/

ml BCFtools/1.15.1-GCC-8.3.0

bcftools stats $BASE/out/9_filt/phylo2_snpsfilt.vcf.gz > $BASE/stats/snps_phylo2_snpsfilt.stats
# SN	0	number of samples:	335
# SN	0	number of records:	45000812
# SN	0	number of no-ALTs:	0
# SN	0	number of SNPs:	45000812
# SN	0	number of MNPs:	0
# SN	0	number of indels:	0
# SN	0	number of others:	0
# SN	0	number of multiallelic sites:	0
# SN	0	number of multiallelic SNP sites:	0

bcftools stats $BASE/out/8_geno/phylo2_all-indel.vcf.gz > $BASE/stats/snps_phylo2_all-indel.stats
# SN	0	number of samples:	335
# SN	0	number of records:	521554932
# SN	0	number of no-ALTs:	445234515
# SN	0	number of SNPs:	74822610
# SN	0	number of MNPs:	0
# SN	0	number of indels:	5989973
# SN	0	number of others:	889415
# SN	0	number of multiallelic sites:	18754526
# SN	0	number of multiallelic SNP sites:	12591365



### ============================================================================
### phyps2 dataset (no Serranus outgroups)
### ============================================================================

BASE=/gss/work/haex1482/phylo2/2_genotyping

grep -Ev "tor|tab|tig" fofn/6_like.fofn > fofn/6_phyps2.fofn   # 329 libraries


### ============================================================================
### 7p. Combine individual GVCF files into cohort GVCF -- phyps2

#!/bin/bash

#SBATCH --job-name=7p_coho
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/7_coho/7p_coho_%j.out
#SBATCH --error=log/7_coho/7p_coho_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
FOFN=$BASE/fofn/6_phyps2.fofn
INPUT=$(cat $FOFN | sed -e 's/^/-V /g')

cd $BASE/out/6_like

gatk --java-options "-Xmx120G" \
    CombineGVCFs \
    $INPUT \
    -O $BASE/out/7_coho/phyps2_cohort_g.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta


### ============================================================================
### 8.1p. Genotype jointly (all sites and SNPs only) -- phyps2

#!/bin/bash

#SBATCH --job-name=8.1p_geno
#SBATCH --partition=carl.p
#SBATCH --array=1-24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/8_geno/8.1p_geno_%A_%a.out
#SBATCH --error=log/8_geno/8.1p_geno_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
SEQ=($(seq -w 1 24))

gatk --java-options "-Xmx120G" \
    GenotypeGVCFs \
    -V $BASE/out/7_coho/phyps2_cohort_g.vcf.gz \
    -L LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]} \
    -O $BASE/out/8_geno/phyps2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --include-non-variant-sites true

# Including indels (executed later)
gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/8_geno/phyps2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -O $BASE/out/8_geno/phyps2_all_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta

gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/8_geno/phyps2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -O $BASE/out/8_geno/phyps2_all-indel_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --select-type-to-exclude INDEL

gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/8_geno/phyps2_intermed_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -O $BASE/out/8_geno/phyps2_rawsnps_LG${SEQ[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --select-type-to-include SNP


### ============================================================================
### 8.2p. Gather LGs and extract variants stats -- phyps2

#!/bin/bash

#SBATCH --job-name=8.2p_geno
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=00-02:00  # D-HH:MM
#SBATCH --output=log/8_geno/8.2p_geno_%j.out
#SBATCH --error=log/8_geno/8.2p_geno_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2/2_genotyping
ALI=$(find $BASE/out/8_geno/phyps2_all_LG*.vcf.gz | awk '{ print "-I", $1 }')
ALL=$(find $BASE/out/8_geno/phyps2_all-indel_LG*.vcf.gz | awk '{ print "-I", $1 }')
RAW=$(find $BASE/out/8_geno/phyps2_rawsnps_LG*.vcf.gz | awk '{ print "-I", $1 }')

# Including indels (executed later)
gatk --java-options "-Xmx120g" \
    GatherVcfs \
    $ALI \
    -O $BASE/out/8_geno/phyps2_all.vcf.gz

tabix -p vcf $BASE/out/8_geno/phyps2_all.vcf.gz

# All sites (except indels)
gatk --java-options "-Xmx120g" \
    GatherVcfs \
    $ALL \
    -O $BASE/out/8_geno/phyps2_all-indel.vcf.gz

tabix -p vcf $BASE/out/8_geno/phyps2_all-indel.vcf.gz

# SNPs only (had to be executed manually)
gatk --java-options "-Xmx120g" \
    GatherVcfs \
    $RAW \
    -O $BASE/out/8_geno/phyps2_rawsnps.vcf.gz

tabix -p vcf $BASE/out/8_geno/phyps2_rawsnps.vcf.gz

# gatk --java-options "-Xmx120G" \
#     VariantsToTable \
#     -V $BASE/out/8_geno/phyps2_rawsnps.vcf.gz \
#     -O $BASE/out/9_filt/phyps2_rawsnps.qual.tsv \
#     -F CHROM -F POS -F MQ -F QD -F FS -F MQRankSum -F ReadPosRankSum \
#     --show-filtered


## -----------------------------------------------------------------------------
## Remove intermediary files (execute manually)
rm $BASE/out/8_geno/phyps2_intermed_*.vcf*
rm $BASE/out/8_geno/phyps2_*LG*


### ============================================================================
### 9. Filter variants -- phyps2

#!/bin/bash

#SBATCH --job-name=9p_filter
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=21-00:00  # D-HH:MM
#SBATCH --output=log/9_filt/9p_filt_%j.out
#SBATCH --error=log/9_filt/9p_filt_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8
ml VCFtools/0.1.16-GCC-8.3.0

BASE=/gss/work/haex1482/phylo2/2_genotyping
NOSAM=$(find $BASE/out/6_like/*.vcf.gz | wc -l)
MISS=$((NOSAM / 10))

gatk --java-options "-Xmx120G" \
    VariantFiltration \
    -V $BASE/out/8_geno/phyps2_rawsnps.vcf.gz \
    -O $BASE/out/9_filt/phyps2_snpsint.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --filter-expression "MQ < 57.5" \
    --filter-name "filter_MQ" \
    --filter-expression "QD < 3.0" \
    --filter-name "filter_QD" \
    --filter-expression "FS > 60.0" \
    --filter-name "filter_FS" \
    --filter-expression "MQRankSum < -0.2 || MQRankSum > 0.2" \
    --filter-name "filter_MQRankSum" \
    --filter-expression "ReadPosRankSum < -2.25 || ReadPosRankSum > 2.25 " \
    --filter-name "filter_ReadPosRankSum"

gatk --java-options "-Xmx120G" \
    SelectVariants \
    -V $BASE/out/9_filt/phyps2_snpsint.vcf.gz \
    -O $BASE/out/9_filt/phyps2_snpsint2.vcf.gz \
    -R $BASE/ref/Hpue_genome_unmasked_01.fasta \
    --exclude-filtered

vcftools \
    --gzvcf $BASE/out/9_filt/phyps2_snpsint2.vcf.gz \
    --max-missing-count $MISS \
    --max-alleles 2 \
    --stdout \
    --recode |
bgzip > $BASE/out/9_filt/phyps2_snpsfilt.vcf.gz

tabix -p vcf $BASE/out/9_filt/phyps2_snpsfilt.vcf.gz


## -----------------------------------------------------------------------------
## Remove intermediary files (execute manually)
rm $BASE/out/9_filt/*snpsint*


### ============================================================================
### Cleanup and stats

mv *.js scripts/

ml BCFtools/1.15.1-GCC-8.3.0

bcftools stats $BASE/out/9_filt/phyps2_snpsfilt.vcf.gz > $BASE/stats/snps_phyps2_snpsfilt.stats
# SN	0	number of samples:	327
# SN	0	number of records:	24794034
# SN	0	number of no-ALTs:	0
# SN	0	number of SNPs:	24794034
# SN	0	number of MNPs:	0
# SN	0	number of indels:	0
# SN	0	number of others:	0
# SN	0	number of multiallelic sites:	0
# SN	0	number of multiallelic SNP sites:	0

bcftools stats $BASE/out/8_geno/phyps2_all-indel.vcf.gz > $BASE/stats/snps_phyps2_all-indel.stats
# SN	0	number of samples:	327
# SN	0	number of records:	544611197
# SN	0	number of no-ALTs:	513064446
# SN	0	number of SNPs:	31145203
# SN	0	number of MNPs:	0
# SN	0	number of indels:	1609881
# SN	0	number of others:	0
# SN	0	number of multiallelic sites:	4601691
# SN	0	number of multiallelic SNP sites:	659900

bcftools stats $BASE/out/8_geno/phyps2_all.vcf.gz > $BASE/out/8_geno/snps_phyps2_all.stats
# SN	0	number of samples:	327
# SN	0	number of records:	551512072
# SN	0	number of no-ALTs:	513064446
# SN	0	number of SNPs:	31145203
# SN	0	number of MNPs:	0
# SN	0	number of indels:	8510756
# SN	0	number of others:	0
# SN	0	number of multiallelic sites:	6385016
# SN	0	number of multiallelic SNP sites:	659900


### ============================================================================
### Soft filtering and stats

# phylo2_all-indel.vcf.gz: 521554932 sites (windows, mtg trees, DILS)
# phyps2_all.vcf.gz: 551512072 sites

# phylo2_snpsfilt.vcf.gz: 335 samples, 45000812 sites
# phyps2_snpsfilt.vcf.gz: 327 samples, 24794034 sites (PCA)

# Filter by minor allele count and physical distance (m2k5)
vcftools \
    --gzvcf phylo2_snpsfilt.vcf.gz \
    --mac 2 \
    --thin 5000 \
    --recode \
    --stdout |
grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phylo2_m2k5.vcf.gz
#> 110436 sites

vcftools \
    --gzvcf phyps2_snpsfilt.vcf.gz \
    --mac 2 \
    --thin 5000 \
    --recode \
    --stdout |
grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phyps_m2k5.vcf.gz
#> 109968 sites

### ----------------------------------------------------------------------------
### Filter effect tests

# Missing data
for i in 0 0.33 0.66 1
do
vcftools \
    --gzvcf /gss/work/haex1482/phylo2/2_genotyping/out/9_filt/phyps2_snpsfilt.vcf.gz \
    --mac 2 \
    --max-missing $i
done

#> 0, 0.33, 0.66: 14612782 out of a possible 24794034 Sites
#>             1: 14028466 out of a possible 24794034 Sites (96% of above)
## Sites already filtered by missing genotypes max 10% of samples

# Minor allele count (phylo2)
for i in 2 4 8
do
vcftools \
    --gzvcf phylo2_snpsfilt.vcf.gz \
    --mac $i
done

#> 2: 35478305 out of a possible 45000812 Sites (79%)
#> 4: 29255378 out of a possible 45000812 Sites (65%)
#> 8: 17666329 out of a possible 45000812 Sites (39%)

# Minor allele count (phyps2)
for i in 2 4 8
do
vcftools \
    --gzvcf phyps2_snpsfilt.vcf.gz \
    --mac $i
done

#> 2: 14612782 out of a possible 24794034 Sites (59%)
#> 4: 10768583 out of a possible 24794034 Sites (43%)
#> 8:  8393043 out of a possible 24794034 Sites (34%)

# Minor allele count and physical distance
for i in 4 6 8 10
do
vcftools \
    --gzvcf phylo2_snpsfilt.vcf.gz \
    --mac $i \
    --thin 5000 \
    --recode \
    --stdout |
grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phylo2_m${i}k5.vcf.gz
done

#> no. of sites ~110 k, independent of mac due to thinning


### ----------------------------------------------------------------------------
### Relabel samples

# Hard-filtered files
bcftools reheader \
  -s relabel.txt \
  --threads 8 \
  phylo2_snpsfilt.vcf.gz |
zcat | grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phylo2e_snpsfilt2.vcf.gz

bcftools reheader \
  -s relabel.txt \
  --threads 8 \
  phyps2_snpsfilt.vcf.gz |
zcat | grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phyps2e_snpsfilt2.vcf.gz

tabix -p vcf phylo2e_snpsfilt2.vcf.gz
tabix -p vcf phyps2e_snpsfilt2.vcf.gz

# Soft-filtered files
bcftools reheader \
  -s relabel.txt \
  --threads 8 \
  phylo2_m2k5.vcf.gz |
zcat | grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phylo2e_m2k5.vcf.gz

bcftools reheader \
  -s relabel.txt \
  --threads 8 \
  phyps2_m2k5.vcf.gz |
zcat | grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phyps2e_m2k5.vcf.gz

tabix -p vcf phylo2e_m2k5.vcf.gz
tabix -p vcf phyps2e_m2k5.vcf.gz

# Relabel id files
cp ids_phylo2.txt ids_phylo2e.txt
cp ids_phyps2.txt ids_phyps2e.txt

while read a b
do
  sed -i "s/$a/$b/" ids_phylo2e.txt
done < relabel.txt

while read a b
do
  sed -i "s/$a/$b/" ids_phyps2e.txt
done < relabel.txt


### ----------------------------------------------------------------------------
### Soft filtering with updated nomenclature

# Physical distance
#> phylo_m2k5.vcf.gz: 110436 sites (IQ-TREE, SVDQuartets)
#> phyps_m2k5.vcf.gz: 109968 sites

# Genetic distance, no missing data, phyps2e (Admixture, D-stats, to do: PCA)
# see dstats_phylo2.sh
#> phyps2e_m2n1.vcf.gz: 14028466 sites
#> phyps2e_m2n1l5.vcf.gz: 935941 sites

# Phased
# see phasing_phylo2.sh
#> phylo2e_phased.vcf.gz (no filters)
#> phylo2_phased_mac2.vcf.gz (minor allele count only)

# Phased, with ancestral allele state, no missing data (tskit)
# see ts_phylo2_aa.sh
#> phylo2e_..._n1.vcf

# Minor allele count only, phyps2e (GWAS)
vcftools \
    --gzvcf phyps2e_snpsfilt.vcf.gz \
    --mac 2 \
    --recode \
    --stdout |
grep -v -e 'ID=Contig' -e '##GATKCommandLine=' |
bgzip > phyps2e_m2.vcf.gz
#> 14612782 sites
