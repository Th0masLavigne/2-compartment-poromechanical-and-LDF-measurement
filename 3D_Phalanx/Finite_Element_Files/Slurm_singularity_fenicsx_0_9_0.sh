#!/bin/bash
### Line above tells to use the bash
### User fill in HERE
### If #SBATCH, not a comment but a directive for the slurm command.
###
### directives:
#
## Required batch arguments
#SBATCH --job-name=ID_63_local_refine
#SBATCH --partition=HighFreq
#SBATCH --ntasks=1
#SBATCH --nodes=1
##SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=32GB
#
## Suggested batch arguments
#SBATCH --mail-type=ALL
#SBATCH --mail-user=thomas.lavigne@ensam.eu
#
## Logging arguments (IMPORTANT)
#SBATCH --output=slurm_%x-%j.out
#SBATCH --error=slurm_%x-%j.err

### Variables Summary
### echo: affiche à l'écran => écrit les informations du job dans le slurm_%x-%j.out
### tout ce qui commence par un "$" est une variable. 
echo ""
echo -e "\033[34m---------------------------------------------------------------------------------------\033[0m"
echo -e "\033[34mVariables Summary: \033[0m"
echo -e "\tWorking Directory = $SLURM_SUBMIT_DIR"
echo -e "\tJob ID = $SLURM_JOB_ID"
echo -e "\tJob Name = $SLURM_JOB_NAME"
echo -e "\tJob Hosts = $SLURM_JOB_NODELIST"
echo -e "\tNumber of Nodes = $SLURM_NNODES"
echo -e "\tNumber of Cores = $SLURM_NTASKS"
echo -e "\tCores per Node = $SLURM_NTASKS_PER_NODE"

### Modules
module purge
echo ""
echo -e "Purge"
module load EasyBuild foss Singularity
echo ""
echo -e "load EasyBuild foss Singularity"
# load matplotlib
module load matplotlib
echo ""
echo -e "load matplotlib"
# 
module load openpyxl/3.0.10
echo ""
echo -e "load openpyxl"
## raccourci pour module load: ml, pour module avail ml av, module list 
# 
## print in the out the versions and available modules 
# module list
## help sur un module:
# module spider SciPy-bundle/2022.05
## pour voir les modules cachés
# module --show hidden avail
# 
cd $SLURM_SUBMIT_DIR
# en sequentiel
# singularity exec /modules/containers/images/dolfinx/dolfinx-0.8.0.sif python3 main_theoretical_cube.py
# en parallele 
mpirun -n $SLURM_NTASKS singularity exec /modules/containers/images/dolfinx/dolfinx-0.9.0.sif python3 main_finger.py
#
### Write any other commands after if needed
# 
# ${SLURM_NTASKS} mettre en argument du fichier python si opti sequentiel/parallèle; par défaut mettre à 1
# https://koor.fr/Python/CodeSamples/SysArgv.wp

### Send Commands
## print the analysis
echo ""
echo -e "Analysis: compute seff $SLURM_JOB_ID"

## seff ${SLURM_JOB_ID}

### EOF
