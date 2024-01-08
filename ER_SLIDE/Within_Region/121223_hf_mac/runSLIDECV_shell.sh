#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --job-name=WR_S_CV
#SBATCH --mail-user=xiaoh@pitt.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mem=100000MB
#SBATCH --cpus-per-task=9
# make sure this is before your main program for it to always run on exit
#trap "echo 'copying files'; rsync -avz * ${SLURM_SUBMIT_DIR}" EXIT

module load gcc/10.2.0
module load r/4.2.0
Rscript SLIDEcv.R  "slide_cv.yaml"