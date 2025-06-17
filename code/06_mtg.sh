### ================================================================================
### Radiation with reproductive isolation in the near-absence of phylogenetic signal
### 06. Mitochondrial tree inference
### By Martin Helmkampf, last edited 2025-06-16
### ================================================================================


### 1. Genotyping the mitogenome

#!/bin/bash

#SBATCH --job-name=1_geno
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --time=0-02:00
#SBATCH --mem-per-cpu=20G
#SBATCH --output=/dev/null
#SBATCH --error=log/sl_geno_%j.err

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=$WORK/phylo2

# phylo2
gatk --java-options "-Xmx20G" \
    GenotypeGVCFs \
    -V $BASE/2_genotyping/out/7_coho/phylo2_cohort_g.vcf.gz \
    -L LG_M \
    -O phylo2_LGM_tmp.vcf.gz \
    -R $BASE/2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
    --include-non-variant-sites true
#> 16856 of 17159 sites (303 sites ungenotyped)

gatk --java-options "-Xmx20G" \
    SelectVariants \
    -V phylo2_LGM_tmp.vcf.gz \
    -O phylo2_LGM_tmp2.vcf.gz \
    -R $BASE/2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
    --select-type-to-exclude INDEL
#> 16502 of 17159 sites (354 indels)

zcat < phylo2_LGM_tmp2.vcf.gz |
    grep -v -e '##contig' -e '##GATKCommandLine=' |
    bgzip > phylo-all_mtg.vcf.gz

tabix -p vcf phylo-all_mtg.vcf.gz

# phyps2
gatk --java-options "-Xmx20G" \
    GenotypeGVCFs \
    -V $BASE/2_genotyping/out/7_coho/phyps2_cohort_g.vcf.gz \
    -L LG_M \
    -O phyps2_LGM_tmp.vcf.gz \
    -R $BASE/2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
    --include-non-variant-sites true
#> 17073 of 17159 sites (86 sites ungenotyped)

gatk --java-options "-Xmx20G" \
    SelectVariants \
    -V phyps2_LGM_tmp.vcf.gz \
    -O phyps2_LGM_tmp2.vcf.gz \
    -R $BASE/2_genotyping/ref/Hpue_genome_unmasked_01.fasta \
    --select-type-to-exclude INDEL
#> 16998 of 17159 sites (75 indels)

zcat < phyps2_LGM_tmp2.vcf.gz |
    grep -v -e '##contig' -e '##GATKCommandLine=' |
    bgzip > phyps2_mtg.vcf.gz

tabix -p vcf phyps2_mtg.vcf.gz

# rm *tmp*


### ============================================================================
### 2. Conversion to fasta

#!/bin/bash

#SBATCH --job-name=2_conv
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --time=0-02:00
#SBATCH --output=/dev/null
#SBATCH --error=log/sl_conv_%j.err

ml hpc-env/8.3
ml VCFtools/0.1.16-GCC-8.3.0

# phylo2
gzip -cd phylo-all_mtg.vcf.gz > phylo-all_mtg.vcf

vcf-to-tab < \
    phylo-all_mtg.vcf |
    sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > \
    phylo-all_mtg.tab

perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl \
    -i phylo-all_mtg.tab > \
    phylo-all_mtg.fas

# phyps2
gzip -cd phyps2_mtg.vcf.gz > phyps2_mtg.vcf

vcf-to-tab < \
    phyps2_mtg.vcf |
    sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > \
    phyps2_mtg.tab

perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl \
    -i phyps2_mtg.tab > \
    phyps2_mtg.fas


# Confirm sequences are of equal length (aligned)
bioawk -c fastx '{ print $name, length($seq) }' phylo-all_mtg.fas   # 16266 aligned sites
bioawk -c fastx '{ print $name, length($seq) }' phyps2_mtg.fas   # 16921 aligned sites

# rm *tab* *.vcf


### ============================================================================
### 3. RAxML combined tree search and bootstrapping analysis under GTR+G model

#!/bin/bash

#SBATCH --job-name=3_raxml
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=5-00:00
#SBATCH --output=log/sl_raxml_%j.out
#SBATCH --error=log/sl_raxml_%j.err

# phylo2
raxml-ng \
    --all \
    --msa phylo-all_mtg.fas \
    --model GTR+G \
    --tree pars{20},rand{20} \
    --bs-trees 200 \
    --threads 12 \
    --worker AUTO \
    --seed 123 \
    --prefix raxml/phylo-all_mtg_GTRG

# phyps2
raxml-ng \
    --all \
    --msa phyps2_mtg.fas \
    --model GTR+G \
    --tree pars{20},rand{20} \
    --bs-trees 200 \
    --threads 12 \
    --worker AUTO \
    --seed 123 \
    --prefix raxml/phyps2_mtg_GTRG
