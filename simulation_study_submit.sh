#!/bin/bash -l

#SBATCH --job-name=simulation_study
#SBATCH --output=simulation_study_%A_%a.out
#SBATCH --error=simulation_study_%A_%a.err

#SBATCH --nodes=1
#SBATCH --ntasks=120
#SBATCH --mem=228g
#SBATCH -t 1:00:00

#SBATCH --account=mfiecas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neher015@umn.edu

cd /home/mfiecas/neher015/human_genetics_final_project

module load R

Rscript simulation_study.R