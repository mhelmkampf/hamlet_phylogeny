#!/bin/bash
# by: Floriane Coulmance: 01/12/2023
# usage:
# 14_gwas.sh -i <PATH> -j <JOB_ID>
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
# Job 0 split genotyping file (hamlets, no outgroups) for 17 species of interest

jobfile0=0_splitvcf.tmp # temp file
cat > $jobfile0 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=0_splitvcf
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/0_splitvcf_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/0_splitvcf_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=06:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load R/4.3.1-foss-2023a
module load BCFtools/1.18-GCC-13.1.0

# Input the bi-allelic genotyping file (hamlets, no outgroups)
INPUT_BI=$BASE_DIR/data/gwas/phyps2e_m2.vcf.gz                   

# Split genotypingg file into 2 for both small & large clades
bcftools view -Oz -S $BASE_DIR/metadata/large_sample_list \${INPUT_BI} > $BASE_DIR/data/large_phyps2e_m2.vcf.gz


EOA



# ------------------------------------------------------------------------------
# Job 1 convert genotyping file to binary plink files for GWAS

jobfile1=1_convert_plink.tmp # temp file
cat > $jobfile1 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=1_convert_plink
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/1_convert_plink_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/1_convert_plink_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=06:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load VCFtools/0.1.16-GCC-13.1.0
module load R/4.3.1-foss-2023a
module load BCFtools/1.18-GCC-13.1.0


# Convert the genotyping file to plink format
vcftools \
      --gzvcf $BASE_DIR/data/large_phyps2e_m2.vcf.gz \
      --plink \
      --out $BASE_DIR/outputs/gxp_clades/large/GxP_plink

# Convert to hapmap / not mandatory to run
plink \
      --file $BASE_DIR/outputs/gxp_clades/large/GxP_plink \
      --recode12 \
      --out $BASE_DIR/outputs/gxp_clades/large/hapmap

# Convert plink genotyping file to binary files to be used in GWAS
plink \
    --noweb \
    --file $BASE_DIR/outputs/gxp_clades/large/GxP_plink \
    --make-bed \
    --out $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary

# Save a copy of the binary .fam file
cp $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary.fam \
   $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary_sauvegarde.fam


EOA



# ------------------------------------------------------------------------------
# Job 2 prepare phenotyping file on species ID (17 species) for GWAS

jobfile2=2_prep.tmp # temp file
cat > $jobfile2 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=2_prep.tmp
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/2_prep_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/2_prep_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=04:00:00


# Input files to merge, binary genotyping file & phenotyping file
fam=$BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary.fam
pheno=$BASE_DIR/metadata/large_pheno.txt


# Sort both input files 
sort -k1 \${fam} > $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary_sorted.fam
sort -k1 \${pheno} > $BASE_DIR/outputs/gxp_clades/large/pheno_sorted

# Merge input files on sampleID
join $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary_sorted.fam $BASE_DIR/outputs/gxp_clades/large/pheno_sorted | \
awk -F " " '{print \$1,\$2,\$3,\$4,\$5,\$8}' \
> $BASE_DIR/outputs/gxp_clades/large/pheno_intermediate

# Format output merged file for use in GWAS
echo -e 'label Within_family_ID ID_father ID_mother Sex spec' \
> $BASE_DIR/outputs/gxp_clades/large/pheno_table.fam && cat $BASE_DIR/outputs/gxp_clades/large/pheno_intermediate \
>> $BASE_DIR/outputs/gxp_clades/large/pheno_table.fam


EOA



# ------------------------------------------------------------------------------
# Job 3 GEMMA univariate analyses (Linear Mixed Models)

jobfile3=3_gemma.tmp # temp file
cat > $jobfile3 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=3_gemma.tmp
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/3_gemma_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/3_gemma_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04:00:00


# Function to read and print to each line
body() {                                                                             
	IFS= read -r header
	printf '%s\n' "\$header"
	"\$@"
}


# Format phenotyping file to match necessary input format & names for GEMMA 
awk '{print \$1, \$2, \$3, \$4, \$5, \$6}' $BASE_DIR/outputs/gxp_clades/large/pheno_table.fam | tail -n +2 > $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary.fam


# GEMMA automatic output folder management
mkdir $BASE_DIR/output/gxp_clades/
mkdir $BASE_DIR/output/gxp_clades/large/


# 1) Relatedness sample matrix
gemma -bfile $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary -gk 1 -o /gxp_clades/large/large

# 2) Linear Mixed Model (lmm)
gemma -bfile $BASE_DIR/outputs/gxp_clades/large/GxP_plink_binary -k output/gxp_clades/large/large.cXX.txt -lmm 4 -o /gxp_clades/large/large.lmm

# 3) Output reformatting (lmm)
sed 's/\\trs\\t/\\tCHROM\\tPOS\\t/g; s/\\([0-2][0-9]\\):/\\1\\t/g' output/gxp_clades/large/large.lmm.assoc.txt | \
      cut -f 2,3,8-10,13-15 | body sort -k1,1 -k2,2n | gzip > $BASE_DIR/outputs/gxp_clades/large/large.lmm.GxP.txt.gz


EOA



# ------------------------------------------------------------------------------
# Job 4 GEMMA averaged results over sliding windows

jobfile4=4_windows.tmp # temp file
cat > $jobfile4 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=4_windows.tmp
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/4_windows_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/4_windows_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=06:00:00


# Input the GEMMA results (lmm)
lmm=$BASE_DIR/outputs/gxp_clades/large/large.lmm.GxP.txt.gz
echo \${lmm}

# First set of parameters, 50kb windows sliding every 5kb 
win5=50000
step5=5000
echo \${win5}
echo \${step5}

# Second set of parameters, 10kb windows sliding every 1kb 
win1=10000
step1=1000
echo \${win1}
echo \${step1}

# Run the average over genome window for each set of parameters
$BASE_DIR/code/sh/gxp_slider.sh \${lm} \${win5} \${step5}
$BASE_DIR/code/sh/gxp_slider.sh \${lm} \${win1} \${step1}
$BASE_DIR/code/sh/gxp_slider.sh \${lmm} \${win5} \${step5}
$BASE_DIR/code/sh/gxp_slider.sh \${lmm} \${win1} \${step1}


EOA



# ------------------------------------------------------------------------------
# Job 5 GWAS plots

jobfile5=5_plots.tmp # temp file
cat > $jobfile5 <<EOA # generate the job file
#!/bin/bash
#SBATCH --job-name=5_plots.tmp
#SBATCH --partition=rosa.p
#SBATCH --output=$BASE_DIR/logs/5_plots_%A_%a.out
#SBATCH --error=$BASE_DIR/logs/5_plots_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=22G
#SBATCH --time=04:00:00


# Load necessary modules for University of Oldenburg cluster (ROSA)
module load hpc-env/13.1
module load R/4.3.1-foss-2023a

# Remove temporary files created during GWAS analyses
rm $BASE_DIR/outputs/gxp_clades/large/\*.tmp

# Run Rscript for GWAS Manhattan plot 
Rscript $BASE_DIR/code/R/gxp_plot.R $BASE_DIR/outputs/gxp_clades/large/ $BASE_DIR/figures/gxp_clades/large/


EOA



# ********** Schedule the job launching ***********
# -------------------------------------------------

if [ "$JID_RES" = "jid1" ] || [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   0_splitvcf      : DONE         **"
else
  jid0=$(sbatch ${jobfile0})
fi


if [ "$JID_RES" = "jid2" ] || [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   1_convert_plink : DONE         **"
elif [ "$JID_RES" = jid1 ]
then
  jid1=$(sbatch ${jobfile1})
else
  jid1=$(sbatch --dependency=afterok:${jid0##* } ${jobfile1})
fi


if [ "$JID_RES" = "jid3" ] || [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   2_prep          : DONE         **"
elif [ "$JID_RES" = "jid2" ]
then
  jid2=$(sbatch ${jobfile2})
else
  jid2=$(sbatch --dependency=afterok:${jid1##* } ${jobfile2})
fi


if [ "$JID_RES" = "jid4" ] || [ "$JID_RES" = "jid5" ];
then
  echo "*****   3_gemma         : DONE         **"
elif [ "$JID_RES" = "jid3" ]
then
  jid3=$(sbatch ${jobfile3})
else
  jid3=$(sbatch --dependency=afterok:${jid2##* } ${jobfile3})
fi


if [ "$JID_RES" = "jid5" ];
then
  echo "*****   4_windows       : DONE         **"
elif [ "$JID_RES" = "jid4" ]
then
  jid4=$(sbatch ${jobfile4})
else
  jid4=$(sbatch --dependency=afterok:${jid3##* } ${jobfile4})
fi


if [ "$JID_RES" = "jid5" ];
then
  jid5=$(sbatch ${jobfile5})
else
  jid5=$(sbatch --dependency=afterok:${jid4##* } ${jobfile5})
fi
