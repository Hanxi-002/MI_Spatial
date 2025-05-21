#!/bin/bash
#SBATCH -t 5-00:00

#SBATCH --job-name=ccr2
# partition (queue) declaration

#SBATCH --mail-user=xiaoh@pitt.edu

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --mem=100g

#SBATCH --cpus-per-task=16

module load gcc/10.2.0
module load r/4.2.0

Rscript run_optimSLIDE.R

