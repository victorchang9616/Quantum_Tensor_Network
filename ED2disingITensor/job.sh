#!/bin/bash 
#SBATCH -J N9-Jpam0.1
#SBATCH -o N9-Jplam0.1.o%j 
#SBATCH --mem 32000
#SBATCH -w, --nodelist=compute-0-3
#SBATCH -N 1
#SBATCH -n 24

module load gcc

./iqdmrg
