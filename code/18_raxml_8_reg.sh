#!/bin/bash
# by: Floriane Coulmance: 12/05/2024
# usage:
# sbatch 18_raxml_8_reg.sh -i <PATH> -j <JOB_ID>
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

while getopts i:j: option
do
case "${option}"
in
i) BASE_DIR=${OPTARG};; # get the base directory path
j) JID_RES=${OPTARG};; # get the jobid from which you want to resume
esac
done



# ********* Create necessary repositories *********
# -------------------------------------------------

# Outputs repo for each step of the pipeline
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
# Job 0 concat the 8 region alignment fasta files into one file

jobfile0=0_concat.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_concat
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/0_concat_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_concat_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=06:00:00


# Use SEQKIT to concatenate the 8 regions fasta alignment files into 1 alignment fasta file
seqkit concat $BASE_DIR/outputs/gxp_clades/large/LG04_1_11555000_11655000.fas $BASE_DIR/outputs/gxp_clades/large/LG08_2_2020000_2080000.fas $BASE_DIR/outputs/gxp_clades/large/LG08_3_23455000_23530000.fas $BASE_DIR/outputs/gxp_clades/large/LG12_4_20135000_20340000.fas $BASE_DIR/outputs/gxp_clades/large/LG12_5_22160000_22255000.fas $BASE_DIR/outputs/gxp_clades/large/LG17_6_22505000_22665000.fas $BASE_DIR/outputs/gxp_clades/large/LG19_7_2375000_2425000.fas $BASE_DIR/outputs/gxp_clades/large/LG23_8_13965000_14050000.fas > $BASE_DIR/outputs/gxp_clades/large/8_regions.fas


EOA



# ------------------------------------------------------------------------------
# Job 1 run raxml on the 8 regions from concatenated fasta file

jobfile1=1_rax.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_rax
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/1_rax_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_rax_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=30G
#SBATCH --time=3-00:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Input the fasta alignment file created at previous step
FAS=$BASE_DIR/outputs/gxp_clades/large/8_regions.fas
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
      --prefix $BASE_DIR/outputs/gxp_clades/large/8_regions


EOA



# ------------------------------------------------------------------------------
# Job 2 plot phylogenetic tree

jobfile2=2_plot.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_plot
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/2_plot_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_plot_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=03:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load R/4.3.1-foss-2023a


# Input the phylogenetic tree support result from previous job
SUP=$BASE_DIR/outputs/gxp_clades/large/8_regions.raxml.support
echo \${SUP}

# Run Rscript to plot phylogenetic tree
Rscript $BASE_DIR/code/R/tree.R $BASE_DIR \${SUP}


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] ||  [ "$JID_RES" = "jid2" ];
then
  echo "*****   0_concat    : DONE         **"
else [ "$JID_RES" = jid0 ]
then
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ];
then
  echo "*****   1_rax      : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid2" ];
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi
