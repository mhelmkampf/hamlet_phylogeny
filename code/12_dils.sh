#!/bin/bash
# by: Floriane Coulmance: 01/12/2023
# usage:
# 12_dils.sh -i <BASE_DIR>
# ------------------------------------------------------------------------------
# <BASE_DIR> is the base directory path, here /hamlet_phylogeny
# All necessary folders and output files will be created by the script
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
esac
done



# ********** Run DILS analysis **********
# -------------------------------------------------

# You need to install the DILS sofware from https://github.com/popgenomics/DILS
# and make sure all the dependency requirements are met https://github.com/dils-popgen/dils/blob/master/manual.pdf  
# Modify the path to DILS_2pop.sh file accordingly to where you have downloaded the DILS software.
# You can run the following line of code directly on your cluster by copy pasting it
# and adding #SBATCH options directly in DILS_2pop.sh, after the shebang

DILS_2pop.sh $BASE_DIR/code/yaml/large_small_ser_2pop.yaml

