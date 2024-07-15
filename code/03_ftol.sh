### ============================================================================
### phylo2
### Phylogenetic reconstruction of Serraninae, including all hamlet species
### ============================================================================

### Preparations

mkdir $WORK/phylo2/3_phylogeny/ftol

base=$WORK/phylo2/3_phylogeny/ftol

cd $base

mkdir 0_prep
mkdir 1_geno
mkdir 2_genes_ftol
mkdir 3_genes_phylo2
mkdir 4_genes_comb
mkdir 5_aln
mkdir 6_iqtree


### ============================================================================
### 0. Previous steps taken (PNAS paper / test run)

# Obtain gene partitions in FToL alignment (0_prep/ftol_partitions.bed < partition.bed)
# Identify genes in reference assembly (0_prep/ftol_coord_Hpue.csv < coord_R24_Hpue_ed.bed)
# FToL per-gene alignments: 2_genes_ftol/*.aln files

# Gblocks / manual alignment editing unnecessary


### ============================================================================
### 1. Re-genotype relevant contigs (LG_M, unplaced)

#!/bin/bash

#SBATCH --job-name=geno
#SBATCH --partition=rosa.p
#SBATCH --array=1-15
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=02-00  # DD-HH
#SBATCH --output=sl_%A_%a_geno.out
#SBATCH --error=sl_%A_%a_geno.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml GATK/4.4.0.0-GCCcore-13.1.0-Java-17

base=$WORK/phylo2/3_phylogeny/ftol
chr=("Contig11544" "Contig11607" "Contig11888" "LG07" "LG12" "LG13" "LG14" "LG16" "LG18" "LG19" "LG20" "LG21" "LG23" "LG24" "LG_M")

gatk --java-options "-Xmx120G" \
    GenotypeGVCFs \
    -V $WORK/phylo2/2_genotyping/out/7_coho/phylo2_cohort_g.vcf.gz \
    -L ${chr[((SLURM_ARRAY_TASK_ID-1))]} \
    -O $base/1_geno/phylo2_allsites_${chr[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $WORK/phylo2/2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
    --include-non-variant-sites true

# Parallel array command block:
gatk --java-options "-Xmx120G" \
    GenotypeGVCFs \
    -V $WORK/phylo2/2_genotyping/out/7_coho/phyps2_cohort_g.vcf.gz \
    -L ${chr[((SLURM_ARRAY_TASK_ID-1))]} \
    -O $base/1_geno/phyps2_allsites_${chr[((SLURM_ARRAY_TASK_ID-1))]}.vcf.gz \
    -R $WORK/phylo2/2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
    --include-non-variant-sites true


### ============================================================================
### 2. Obtain FToL per-gene alignments

# copy from PNAS project > 2_genes_ftol/*.aln files


### ============================================================================
### 3.1. Select samples and extract sites

$base/0_prep/samples_ingroup.ids

## Max coverage unless noted otherwise
# abe PL17_88abebel (pnas: 20864abehon)
# aff 19174affgun
# atl 54786atlliz
# cas 54649casliz (mt: Caribbean)
# chl PL17_38chlpri
# eco 62576ecoarc
# flo 62558floarc (pnas: PL17_160floflk, mt: gulf)
# gem 62570gemarc
# gum 23301gumboc (pnas: 20642gumhon)
# gut 19076gutbar (2nd highest cov, 1st is outlier)
# ind PL17_64indpri (pnas: 18238indbel)
# lib HypoHaiti1libhai
# may PL17_122maybel (pnas: same)
# nig 18906nigboc (highest cov of major pop; pnas: 18906nigpan)
# pro 20650prohon
# pue 19104puegun (pnas: 18434puepan)
# ran FL0880ransan (pnas: 20613ranhon) 
# uni PL17_136uniflk (pnas: 18448unipan)

$base/0_prep/samples_outgroup.ids

# tab 20480tabhon
# tig PL17_21tigboc
# tor Bocas16.3torboc

# v1: also includes tan PL17_72tanpri


#!/bin/bash

#SBATCH --job-name=extract
#SBATCH --partition=rosa.p
#SBATCH --array=1-23
#SBATCH --nodes=1
#SBATCH --time=00-01  # DD-HH
#SBATCH --output=sl_%A_%a_extract.out
#SBATCH --error=sl_%A_%a_extract.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml VCFtools/0.1.16-GCC-13.1.0

base=$WORK/phylo2/3_phylogeny/ftol
lines=$(head $base/0_prep/ftol_coord_Hpue.csv -n $SLURM_ARRAY_TASK_ID | tail -n 1)
read chr start end gene dir <<< $lines

vcftools --gzvcf $base/1_geno/phyps2_allsites_${chr}.vcf.gz \
    --keep $base/0_prep/samples_ingroup.ids \
    --chr $chr \
    --from-bp $start \
    --to-bp $end \
    --recode \
    --stdout |
grep -v "##" > $base/3_genes_phylo2/${gene}_in.vcf

vcftools --gzvcf $base/1_geno/phylo2_allsites_${chr}.vcf.gz \
    --keep $base/0_prep/samples_outgroup.ids \
    --chr $chr \
    --from-bp $start \
    --to-bp $end \
    --recode \
    --stdout |
grep -v "##" > $base/3_genes_phylo2/${gene}_out.vcf


### ============================================================================
### 3.2. Convert to fasta, join and reorient (*** command line ***)

ml VCFtools/0.1.16-GCC-13.1.0

base=$WORK/phylo2/3_phylogeny/ftol
genes=`awk '{ print $4 }' $base/0_prep/ftol_coord_Hpue.csv`
rev=`awk '$5 == "-" { print $4 }' $base/0_prep/ftol_coord_Hpue.csv`

for file in $base/3_genes_phylo2/*.vcf
do
    vcf-to-tab < $file > ${file%.vcf}.tab
    perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl -i ${file%.vcf}.tab \
    > ${file%.vcf}.fas
    rm ${file%.vcf}.tab*
done

for i in $genes
do
    cat $base/3_genes_phylo2/${i}_in.fas $base/3_genes_phylo2/${i}_out.fas \
    > $base/3_genes_phylo2/${i}_phylo2.fas
done

for i in $rev
do
    mv $base/3_genes_phylo2/${i}_phylo2.fas $base/3_genes_phylo2/${i}_phylo2.fas.org
    seqtk seq -r $base/3_genes_phylo2/${i}_phylo2.fas.org \
    > $base/3_genes_phylo2/${i}_phylo2.fas
done


### ============================================================================
### 4. Unalign FToL data and combine with phylo2 data (*** command line ***)

base=$WORK/phylo2/3_phylogeny/ftol
genes=`awk '{ print $4 }' $base/0_prep/ftol_coord_Hpue.csv`

$base/0_prep/taxa_to_exclude.ids
# Hypoplectrus_aberrans
# Hypoplectrus_chlorurus
# Hypoplectrus_gemma
# Hypoplectrus_gummigutta
# Hypoplectrus_guttavarius
# Hypoplectrus_indigo
# Hypoplectrus_nigricans
# Hypoplectrus_puella
# Hypoplectrus_unicolor
# Serranus_tortugarum
# Serranus_tabacarius
# Serranus_tigrinus

cd $base

for i in $genes
do
    sed -e 's/N//g' -e 's/> />/g' -e '/^$/d' $base/2_genes_ftol/${i}_ftol.aln | 
    awk 'BEGIN { while ((getline < "0_prep/taxa_to_exclude.ids") >0) l[">"$1]=1 } /^>/ { f=!l[$1] }f' \
    > $base/2_genes_ftol/${i}_ftol.fas
done


for i in $genes
do
    cat $base/3_genes_phylo2/${i}_phylo2.fas $base/2_genes_ftol/${i}_ftol.fas \
    > $base/4_genes_comb/${i}_comb.fas
done


### ============================================================================
### 5. Realign and edit (*** command line, ca. 4 min ***)

ml MAFFT/7.505-GCC-13.1.0-with-extensions

base=$WORK/phylo2/3_phylogeny/ftol
genes=`awk '{ print $4 }' $base/0_prep/ftol_coord_Hpue.csv`

cd $base

for i in $genes
do
    mafft --maxiterate 1000 --globalpair --adjustdirectionaccurately $base/4_genes_comb/${i}_comb.fas \
    > $base/5_aln/${i}_comb.aln
done

# first sample must have sequence data
# visual inspection: overall well aligned between ftol-phylo2 species pairs, but some inconsistent indels


## Gblocks not necessary, only applicable to 4 genes
# for i in 12s 16s coi cytb
# do
#    Gblocks $base/5_aln/v1/${i}_comb.aln -t=D -b2=0 -b3=4 -b4=5 -b5=a -e=.gb -v=100
#    mv $base/5_aln/v1/${i}_comb.aln.gb $base/5_aln/v2/${i}_comb.aln.gb
#    rm $base/5_aln/v2/${i}_comb.aln
# done

# 12s:  97%
# 16s:  75%
# coi:  97%
# cytb: 69%


### ============================================================================
### 6. Phylogenetic inference (IQ-TREE)

#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH --partition=rosa.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01-00  # DD-HH
#SBATCH --output=sl_%j_iqtree.out
#SBATCH --error=sl_%j_iqtree.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@uni-oldenburg.de

ml IQ-TREE/2.2.2.7-gompi-2023a

base=$WORK/phylo2/3_phylogeny/ftol

cd $base/6_iqtree

iqtree2 \
    -p $base/5_aln \
    -b 200 \
    -allnni \
    -T AUTO \
    --prefix ftol_v2


# grep 'parsimony-informative' *.log | awk -F " " '{ sum += $1 } END { print sum }'
#> 1544

# grep 'singleton sites' *.log | awk -F " " '{ sum += $3 } END { print sum }'
#> 1235

# grep 'constant sites' *.log | awk -F " " '{ sum += $6 } END { print sum }'
#> 15038

# partition model based on gene-specific model tests


## ----------------------------------------------------------------------------
## x. New versioning

# 6_iqtree/allftol_included ("all"): no tan, ftol samples included, no alignment editing
# v2: no tan, redundant ftol samples excluded, no alignment editing


## ----------------------------------------------------------------------------
## x. Old versioning (v1)

# v1: no ftol samples excluded, no alignment editing

# v2: as v1, but Gblocks editing of 12s 16s coi cytb
# cd $base/6_iqtree/v2 && iqtree2 -p $base/5_aln/v2 -B 1000 -T AUTO --prefix r23_v2

# v3: as v1, but coi cytb removed manually from phylo2 Serranus
# cd $base/6_iqtree/v3 && iqtree2 -p $base/5_aln/v3 -B 1000 -T AUTO --prefix r23_v3

# v4: removed ftol samples in taxa_to_exclude_v4.ids
# cd $base/6_iqtree/v4 && iqtree2 -p $base/5_aln/v4 -B 1000 -T AUTO --prefix r23_v4

# v5: as v4, but with floflk instead of floarc?