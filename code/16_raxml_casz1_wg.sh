#!/bin/bash
# by: Floriane Coulmance: 17/03/2023
# usage:
# sbatch 16_raxml_casz1_wg.sh -i <PATH> -j <JOB_ID> 
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

# Repo for job logs
mkdir $BASE_DIR/logs/

# Repo for figures
mkdir $BASE_DIR/figures/
mkdir $BASE_DIR/figures/casz1_wholegene_phylo/

# Outputs repo for each step of the pipeline
mkdir $BASE_DIR/outputs/
mkdir $BASE_DIR/outputs/casz1_wholegene_phylo/
mkdir $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/



# ********* Jobs creation *************************
# -------------------------------------------------

# ------------------------------------------------------------------------------
# Job 0 subset genotyping all indels file (hamlets, and outgroups)
# for casz1 gene entirely

jobfile0=0_gen.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_gen
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/0_gen_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_gen_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=10:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Input the all indel genotyping file (hamlets, and outgroups)
VCF=$BASE_DIR/data/trees/phylo2_all-indel.vcf.gz

# Extract Linkage Group (LG), start and end positions corresponding to the gene of interest (casz1)
ANNO=LG12
echo \${ANNO}

START=\$(zless ~/data/annotations/HP.annotation.named.\${ANNO}.gff.gz | grep -w gene | grep -i casz1 | awk '{print \$4}') # get the start position
END=\$(zless ~/data/annotations/HP.annotation.named.\${ANNO}.gff.gz | grep -w gene | grep -i casz1  | awk '{print \$5}') # get the end position
echo \${START}
echo \${END}

# Subset all indel genotyping file based on gene of interest LG, start and end positions
vcftools --gzvcf \${VCF} --chr \${ANNO} --from-bp \${START} --to-bp \${END} --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.vcf.gz

# Create coordinate file to the subsetted gene genotyping file
tabix -p vcf $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.vcf.gz

EOA



# ------------------------------------------------------------------------------
# Job 1 convert genotyping file to fasta alignment

jobfile1=1_tab.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_tab
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/1_tab_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_tab_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=1-00:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Unzip the vcf file for reading
gunzip $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.vcf.gz

# Convert vcf to fasta file in 2 steps, 2nd step uses custom script from https://github.com/JinfengChen/vcf-tab-to-fasta (you need to download this to your base directory)
vcf-to-tab < $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.vcf | sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.tab
perl $BASE_DIR/vcf_tab_to_fasta_alignment.pl -i $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.tab > $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.fas


EOA



# ------------------------------------------------------------------------------
# Job 2 run raxml phylogenetic analysis

jobfile2=2_rax.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_rax
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/2_rax_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_rax_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
##SBATCH --mem-per-cpu=30G
#SBATCH --time=4-00:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Input the corresponding region fasta alignment file created at previous step
FAS=$BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.fas
echo \${FAS}

# Set the phylogenetic model to be considered for phylogenetic reconstruction
MODEL="GTR+G"
echo \${MODEL}

# Reconstruct phylogeny
    ~/apps/raxml-ng/bin/raxml-ng --all \
      --msa \${FAS} \
      --model \${MODEL}  \
      --tree pars{20},rand{20} \
      --bs-trees 100 \
      --threads 24 \
      --seed 123 \
      --prefix $BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene


EOA



# ------------------------------------------------------------------------------
# Job 3 plot phylogenetic tree

jobfile3=3_plot.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_plot
#SBATCH --partition=rosa.p
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


# Input the phylogenetic tree support result from previous job
SUP=$BASE_DIR/outputs/casz1_wholegene_phylo/vcf_casz1_wholegene/casz1_wholegene.raxml.support
echo \${SUP}

# Run Rscript to plot phylogenetic tree
Rscript $BASE_DIR/figures/Fig3c_S10-19_trees.R $BASE_DIR \${SUP}


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
