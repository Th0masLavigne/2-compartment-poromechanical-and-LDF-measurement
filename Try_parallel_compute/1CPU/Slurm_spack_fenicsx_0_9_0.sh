#!/bin/bash

#SBATCH --job-name=1CPU
#SBATCH --partition=CPU_Compute
#SBATCH --ntasks=1
#SBATCH --nodes=1
##SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=16GB
##SBATCH --mem=0
#
## Suggested batch arguments
#SBATCH --mail-type=ALL
#SBATCH --mail-user=thomas.lavigne@ensam.eu
#
## Logging arguments (IMPORTANT)
#SBATCH --output=slurm_%x-%j.out
#SBATCH --error=slurm_%x-%j.err


. /modules/spack/v0.23/share/spack/setup-env.sh
spack env activate tlavigne

srun python main_finger.py

