#!/bin/bash
#SBATCH --job-name=fig1a_R
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1
#SBATCH --time=7-00:00:00
#SBATCH --output=./fig1a_out.%j
#SBATCH --error=./fig1a_err.%j

cd /home/nieduszynski/rosemary/MammalianNanopore/Figures

echo "PATH = " $PATH

PATH=/home/nieduszynski/rosemary/R/x86_64-pc-linux-gnu-library/3.5:$PATH
PATH=/software/TOOLS/R/R-3.5.1/library:$PATH

echo "PATH = " $PATH

# This may vary per HPC system. At USC's hpc system
# we use: source /usr/usc/R/default/setup.sh
#module load R

Rscript test.R /data/nieduszynski_2/workspace/rose/2018_09_10_RW_ONT_RPE_24h_2h/all /data/nieduszynski_2/workspace/rose/2018_09_10_RW_ONT_RPE_24h_2h/all/figure1a_all barcode06,barcode07,barcode05,barcode09 barcode06,barcode10,barcode11,barcode12
