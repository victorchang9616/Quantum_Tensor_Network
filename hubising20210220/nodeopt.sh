#!/bin/bash 
#SBATCH -J TLL
#SBATCH --output=LLN64Szr0.1_%A_%a.txt
##SBATCH --nodelist=compute-0-4
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem 200000
#SBATCH --array=1-30
##SBATCH -t 13-00:00:00
##SBATCH --mail-type=END
##SBATCH --mail-user=schang21@Central.uh.edu
module load intel
module load gcc/6.4.0

./exthubbard inputfile_optlat
