#!/bin/bash 
#SBATCH -J TLL
#SBATCH --output=3x3hz1_%A.txt
##SBATCH --output=2x2nohz%A_%a.txt
##SBATCH --nodelist=compute-0-4
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem 220000
##SBATCH --array=1-100
##SBATCH -t 13-00:00:00
##SBATCH --mail-type=END
##SBATCH --mail-user=schang21@Central.uh.edu
module load intel
module load gcc/6.4.0

./iqdmrg input_ed
