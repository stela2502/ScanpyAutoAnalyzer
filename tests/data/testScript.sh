#! /bin/bash
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -J SingSing
#SBATCH -o SingSing.%j.out
#SBATCH -e SingSing.%j.err
module purge
module load Singularity/default SingSingCell/1.1 
exit 0

