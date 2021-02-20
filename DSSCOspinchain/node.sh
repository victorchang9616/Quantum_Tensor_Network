#!/bin/bash 
#SBATCH -J biqAF
#SBATCH -o 3hes.o%j
##SBATCH --nodelist=compute-0-5
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem 200000
##SBATCH -t 13-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=schang21@Central.uh.edu
module load intel
module load gcc/6.4.0

./dmrg input 
