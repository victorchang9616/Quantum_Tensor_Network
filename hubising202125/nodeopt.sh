#!/bin/bash 
#SBATCH -J TLL
#SBATCH --output=TLLN128nr0.01.%A_%a.txt
##SBATCH --nodelist=compute-0-4
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem 60000
#SBATCH --array=1-10
##SBATCH -t 13-00:00:00
##SBATCH --mail-type=END
##SBATCH --mail-user=schang21@Central.uh.edu
module load intel
module load gcc/6.4.0

./exthubbard inputfile_optlat
