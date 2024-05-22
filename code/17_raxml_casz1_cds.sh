#!/bin/bash
# by: Floriane Coulmance: 17/03/2023
# usage:
# sbatch 17_raxml_casz1_cds.sh -i <PATH> -j <JOB_ID> 
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
mkdir $BASE_DIR/figures/casz1_cds_phylo/

# Outputs repo for each step of the pipeline
mkdir $BASE_DIR/outputs/
mkdir $BASE_DIR/outputs/casz1_cds_phylo/
mkdir $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/



# ********* Jobs creation *************************
# -------------------------------------------------

# ------------------------------------------------------------------------------
# Job 0 subset genotyping all indels file (hamlets, and outgroups)
# for the casz1 gene 23 CDS sequences (exons)

jobfile0=0_prep.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_prep
#SBATCH --partition=rosa.p
#SBATCH --array=1-23
#SBATCH --output=$BASE_DIR/logs/0_prep_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_prep_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=10:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
# module load R/4.3.1-foss-2023a
module load BCFtools/1.18-GCC-13.1.0

# Input the all indel genotyping file (hamlets, and outgroups)
VCF=$BASE_DIR/data/trees/phylo2_all-indel.vcf.gz


# Run this job in parallel for all the 23 casz1 exons
# Subset the all indel genotyping file into 1 file for each of the 23 casz1 exons

if [ "\${SLURM_ARRAY_TASK_ID}" = "1" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20186151 --to-bp 20186375 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon1_20186151_20186375.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "2" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20187957 --to-bp 20188013 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon2_20187957_20188013.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "3" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20216287 --to-bp 20216358 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon3_20216287_20216358.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "4" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20258206 --to-bp 20258253 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon4_20258206_20258253.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "5" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20274561 --to-bp 20274958 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon5_20274561_20274958.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "6" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20316651 --to-bp 20317196 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon6_20316651_20317196.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "7" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20322900 --to-bp 20323800 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon7_20322900_20323800.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "8" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20325613 --to-bp 20325678 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon8_20325613_20325678.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "9" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20328266 --to-bp 20328356 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon9_20328266_20328356.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "10" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20330744 --to-bp 20330908 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon10_20330744_20330908.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "11" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20332529 --to-bp 20332701 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon11_20332529_20332701.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "12" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20333323 --to-bp 20334191 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon12_20333323_20334191.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "13" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20337811 --to-bp 20337943 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon13_20337811_20337943.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "14" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20338356 --to-bp 20338419 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon14_20338356_20338419.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "15" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20339401 --to-bp 20339549 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon15_20339401_20339549.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "16" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20339836 --to-bp 20339958 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon16_20339836_20339958.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "17" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20340064 --to-bp 20340384 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon17_20340064_20340384.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "18" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20340481 --to-bp 20340679 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon18_20340481_20340679.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "19" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20341505 --to-bp 20341676 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon19_20341505_20341676.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "20" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20343248 --to-bp 20343396 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon20_20343248_20343396.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "21" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20344386 --to-bp 20344599 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon21_20344386_20344599.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "22" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20345420 --to-bp 20345785 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon22_20345420_20345785.vcf.gz

elif [ "\${SLURM_ARRAY_TASK_ID}" = "23" ]
then
  vcftools --gzvcf \${VCF} --chr LG12 --from-bp 20346547 --to-bp 20347811 --recode --stdout | bgzip > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon23_20346547_20347811.vcf.gz

else
  echo "please verify #SBATCH-array param"
fi


EOA



# ------------------------------------------------------------------------------
# Job 1 concatenate the 23 exons subsetted all indel genotyping files

jobfile1=1_gen.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_gen
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/1_gen_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_gen_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=10:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Use BCFTOOlS to concatenate the 23 files
bcftools concat $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon1_20186151_20186375.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon2_20187957_20188013.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon3_20216287_20216358.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon4_20258206_20258253.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon5_20274561_20274958.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon6_20316651_20317196.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon7_20322900_20323800.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon8_20325613_20325678.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon9_20328266_20328356.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon10_20330744_20330908.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon11_20332529_20332701.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon12_20333323_20334191.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon13_20337811_20337943.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon14_20338356_20338419.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon15_20339401_20339549.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon16_20339836_20339958.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon17_20340064_20340384.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon18_20340481_20340679.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon19_20341505_20341676.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon20_20343248_20343396.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon21_20344386_20344599.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon22_20345420_20345785.vcf.gz $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_exon23_20346547_20347811.vcf.gz -Oz -o $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.vcf.gz

# Create coordinate file to concatenated 23 exons file
tabix -p vcf $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.vcf.gz


EOA



# ------------------------------------------------------------------------------
# Job 2 convert genotyping file to fasta alignment

jobfile2=2_tab.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_tab
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/2_tab_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_tab_%A_%a.err
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
gunzip $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.vcf.gz

# Convert vcf to fasta file in 2 steps, 2nd step uses custom script from https://github.com/JinfengChen/vcf-tab-to-fasta (you need to download this to your base directory)
vcf-to-tab < $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.vcf | sed -e 's/\.\/\./N\/N/g' -e 's/[ACGTN\*]\/\*/N\/N/g' -e 's/\*\/[ACGTN\*]/N\/N/g' > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.tab
perl $BASE_DIR/vcf_tab_to_fasta_alignment.pl -i $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.tab > $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.fas


EOA



# ------------------------------------------------------------------------------
# Job 3 run raxml phylogenetic analysis

jobfile3=3_rax.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_rax
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/3_rax_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_rax_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
##SBATCH --mem-per-cpu=30G
#SBATCH --time=4-00:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load BCFtools/1.18-GCC-13.1.0

# Input the corresponding fasta alignment file created at previous step
FAS=$BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.fas
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
      --prefix $BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons


EOA



# ------------------------------------------------------------------------------
# Job 4 plot phylogenetic tree

jobfile4=4_plot.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_plot
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/4_plot_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_plot_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90G
#SBATCH --time=03:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load R/4.3.1-foss-2023a


# Input the phylogenetic tree support result from previous job
SUP=$BASE_DIR/outputs/casz1_cds_phylo/vcf_casz1_cds/casz1_23exons.raxml.support
echo \${SUP}

# Run Rscript to plot phylogenetic tree
Rscript $BASE_DIR/figures/Fig3c_S10-19.R $BASE_DIR \${SUP}


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] ||  [ "$JID_RES" = "jid3" ] ||  [ "$JID_RES" = "jid4" ];
then
  echo "*****   0_prep      : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] ||  [ "$JID_RES" = "jid3" ] ||  [ "$JID_RES" = "jid4" ];
then
  echo "*****   1_gen      : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] ||  [ "$JID_RES" = "jid4" ];
then
  echo "*****   2_tab      : DONE         **"
elif [ "$JID_RES" = jid2 ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ];
then
  echo "*****   3_rax      : DONE         **"
elif [ "$JID_RES" = jid3 ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid4" ];
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi
