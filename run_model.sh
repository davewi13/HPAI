#!/bin/sh

#SBATCH --job-name=HPAI_RJMCMC
#SBATCH --partition=medium
#SBATCH --array=1-4
#SBATCH --mem=1G

./run_model.out $SLURM_ARRAY_TASK_ID