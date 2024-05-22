#!/bin/bash
# by: Floriane Coulmance: 01/12/2023
# usage:
# 15_raxml_reg.sh -i <PATH> -j <JOB_ID>
# ------------------------------------------------------------------------------
# <PATH> corresponds to the path to the base directory, all outputs and necessary
# folders will be created by the script
# <JOB_ID> corresponds to ID of the job you want to run this script from within
# this pipeline/file
# ------------------------------------------------------------------------------

#SBATCH --partition=rosa.p
#SBATCH --nodes=1
#SBATCH --time=01:00:00



# ********** Allow to enter bash options **********
# -------------------------------------------------

while getopts i:j:k: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) JID_RES=${OPTARG};; # get the jobid from which you want to resume
esac
done



# ********* Create necessary repositories *********
# -------------------------------------------------

# Repo for gxp outputs
mkdir $BASE_DIR/outputs/
mkdir $BASE_DIR/outputs/gxp_clades/
mkdir $BASE_DIR/outputs/gxp_clades/large/

# Repo for the figures
mkdir $BASE_DIR/figures/
mkdir $BASE_DIR/figures/gxp_clades/
mkdir $BASE_DIR/figures/gxp_clades/large/



# ********* Jobs creation *************************
# -------------------------------------------------

# ------------------------------------------------------------------------------
# Job 0 subset genotyping all indels file (hamlets, and outgroups)
# for each region of interest identified by GWAS (8)

jobfile0=0_gen.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_gen
#SBATCH --partition=rosa.p
#SBATCH --array=0-7
#SBATCH --output=$BASE_DIR/logs/0_gen_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_gen_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=5:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Input the all indel genotyping file (hamlets, and outgroups)
VCF=$BASE_DIR/data/trees/phylo2_all-indel.vcf.gz
echo \${VCF}

# Input file with list and genomic positions of regions of interest from GWAS
INPUT=$BASE_DIR/outputs/gxp_clades/large/gxp_large_gemma_lmm_regions.txt
echo \${INPUT}

# Extract region ID column (8 lines)
REG=(\$(awk '{print \$1}' \${INPUT}))
echo \${REG[@]}

# Run this job in parallel for each region ID (8 jobs)
GEN=\${REG[\${SLURM_ARRAY_TASK_ID}]}
echo \${GEN}

# Extract line from region table corresponding to current region
LINE=(\$(grep "\${GEN}" \${INPUT}))
echo \${LINE[@]}

# Extract Linkage Group (LG), start and end positions corresponding to region considered in this job
LG="\${GEN:0:4}" 
START=\${LINE[2]}
END=\${LINE[3]}

echo \${LG}
echo \${START}
echo \${END}

# Subset all indel genotyping file based on LG, start and end position
vcftools --gzvcf \${VCF} --chr \${LG} --from-bp \${START} --to-bp \${END} --recode --stdout | bgzip > $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.vcf.gz

# Create coordinate file to the subsetted region output file
tabix -p vcf $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.vcf.gz


EOA



# ------------------------------------------------------------------------------
# Job 1 convert each region all indel genotyping file to fasta alignment

jobfile1=1_tab.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_tab
#SBATCH --partition=rosa.p
#SBATCH --array=0-7
#SBATCH --output=$BASE_DIR/logs/1_tab_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_tab_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=5:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0


# Input file with list and genomic positions of regions of interest from GWAS
INPUT=$BASE_DIR/outputs/gxp_clades/large/gxp_large_gemma_lmm_regions.txt
echo \${INPUT}

# Extract region ID column (8 lines)
REG=(\$(awk '{print \$1}' \${INPUT}))
echo \${REG[@]}

# Run this job in parallel for each region ID (8 jobs)
GEN=\${REG[\${SLURM_ARRAY_TASK_ID}]}
echo \${GEN}

# Extract line from region table corresponding to current region
LINE=(\$(grep "\${GEN}" \${INPUT}))
echo \${LINE[@]}

# Extract Linkage Group (LG), start and end positions corresponding to region considered in this job
LG="\${GEN:0:4}" 
START=\${LINE[2]}
END=\${LINE[3]}

echo \${LG}
echo \${START}
echo \${END}

# Input the corresponding region genotyping file created at the previous step
VCF=$BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.vcf.gz
echo \${VCF}

# Unzip the vcf file for reading
gunzip \${VCF}

# Convert vcf to fasta file in 2 steps, 2nd step uses custom script from https://github.com/JinfengChen/vcf-tab-to-fasta (you need to download this to your base directory)
vcf-to-tab < $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.vcf | sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.tab
perl $BASE_DIR/vcf_tab_to_fasta_alignment.pl -i $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.tab > $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.fas


EOA



# ------------------------------------------------------------------------------
# Job 2 run raxml phylogenetic analysis on each region fasta alignement

jobfile2=2_rax.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_rax
#SBATCH --partition=rosa.p
#SBATCH --array=0-7
#SBATCH --output=$BASE_DIR/logs/2_rax_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_rax_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=30G
#SBATCH --time=3-00:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0


# Input file with list and genomic positions of regions of interest from GWAS
INPUT=$BASE_DIR/outputs/gxp_clades/large/gxp_large_gemma_lmm_regions.txt
echo \${INPUT}

# Extract region ID column (8 lines)
REG=(\$(awk '{print \$1}' \${INPUT}))
echo \${REG[@]}

# Run this job in parallel for each region ID (8 jobs)
GEN=\${REG[\${SLURM_ARRAY_TASK_ID}]}
echo \${GEN}

# Extract line from region table corresponding to current region
LINE=(\$(grep "\${GEN}" \${INPUT}))
echo \${LINE[@]}

# Extract Linkage Group (LG), start and end positions corresponding to region considered in this job
LG="\${GEN:0:4}" 
START=\${LINE[2]}
END=\${LINE[3]}

echo \${LG}
echo \${START}
echo \${END}

# Set the phylogenetic model to be considered for phylogenetic reconstruction
MODEL="GTR+G"
echo \${MODEL}

# Input the corresponding region fasta alignment file created at previous step
FAS=$BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.fas
echo \${FAS}

# Reconstruct phylogeny
    ~/apps/raxml-ng/bin/raxml-ng --all \
      --msa \${FAS} \
      --model \${MODEL}  \
      --tree pars{20},rand{20} \
      --bs-trees 100 \
      --threads 24 \
      --seed 123 \
      --prefix $BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}


EOA



# ------------------------------------------------------------------------------
# Job 3 plot phylogenetic tree

jobfile3=3_plot.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_plot
#SBATCH --partition=rosa.p
#SBATCH --array=0-7
#SBATCH --output=$BASE_DIR/logs/3_plot_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_plot_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=03:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load R/4.3.1-foss-2023a


# Input file with list and genomic positions of regions of interest from GWAS
INPUT=$BASE_DIR/outputs/gxp_clades/large/gxp_large_gemma_lmm_regions.txt
echo \${INPUT}

# Extract region ID column (8 lines)
REG=(\$(awk '{print \$1}' \${INPUT}))
echo \${REG[@]}

# Run this job in parallel for each region ID (8 jobs)
GEN=\${REG[\${SLURM_ARRAY_TASK_ID}]}
echo \${GEN}

# Extract line from region table corresponding to current region
LINE=(\$(grep "\${GEN}" \${INPUT}))
echo \${LINE[@]}

# Extract Linkage Group (LG), start and end positions corresponding to region considered in this job
LG="\${GEN:0:4}" 
START=\${LINE[2]}
END=\${LINE[3]}

echo \${LG}
echo \${START}
echo \${END}

# Input the phylogenetic tree support result from previous job
SUP=$BASE_DIR/outputs/gxp_clades/large/\${GEN}_\${START}_\${END}.raxml.support
echo \${SUP}

# Run Rscript to plot phylogenetic tree
Rscript $BASE_DIR/figures/Fig3c_S10-19.R $BASE_DIR \${SUP}


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ];
then
  echo "*****   0_gen      : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ];
then
  echo "*****   1_tab      : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ];
then
  echo "*****   2_rax      : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid3" ];
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi
