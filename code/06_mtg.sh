### ============================================================================
### phylo2
### Mitochondrial genome reconstruction and phylogeny
### ============================================================================

### 1. Genotyping the mitogenome

#!/bin/bash

#SBATCH --job-name=1_geno
#SBATCH --partition=carl.p
#SBATCH --nodes=1
#SBATCH --time=0-02:00
#SBATCH --mem-per-cpu=20G
#SBATCH --output=/dev/null
#SBATCH --error=log/sl_geno_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8

BASE=/gss/work/haex1482/phylo2

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
    bgzip > phylo2_mtg.vcf.gz

tabix -p vcf phylo2_mtg.vcf.gz

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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

ml hpc-env/8.3
ml VCFtools/0.1.16-GCC-8.3.0

# phylo2
gzip -cd phylo2_mtg.vcf.gz > phylo2_mtg.vcf

vcf-to-tab < \
    phylo2_mtg.vcf |
    sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > \
    phylo2_mtg.tab

perl ~/apps/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl \
    -i phylo2_mtg.tab > \
    phylo2_mtg.fas

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
bioawk -c fastx '{ print $name, length($seq) }' phylo2_mtg.fas   # 16266 aligned sites
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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=martin.helmkampf@leibniz-zmt.de

# phylo2
raxml-ng \
    --all \
    --msa phylo2_mtg.fas \
    --model GTR+G \
    --tree pars{20},rand{20} \
    --bs-trees 200 \
    --threads 12 \
    --worker AUTO \
    --seed 123 \
    --prefix raxml/phylo2_mtg_GTRG

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


## ----------------------------------------------------------------------------
## x. Identity / similarity

# Proxy samples
# 18172puebel: Caribbean top haplotype group
# 18428nigboc: Caribbean top haplotype group
# 62558floarc: Gulf haplotype group
# Bocas16.3torboc: Se. outgroup

# Pairwise ident / sim calculation at https://www.bioinformatics.org/sms2/ident_sim.html

#> Within Caribbean top haplotype group: ~40 substitutions (99.7 % ident)
#> Caribbean vs Gulf haplotype groups:  ~650 substitutions (96.0 % ident)
#> Caribbean vs outgroup haplotypes:   ~4300 substitutions (73.0 % ident)
#? Role of Ns?
