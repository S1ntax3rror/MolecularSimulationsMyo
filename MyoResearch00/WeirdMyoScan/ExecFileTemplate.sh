#!/bin/bash
header to replace
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=1000

#--------------------------------------
# Modules
#--------------------------------------

my_orca=/data/kaeserj/GaussianOpt/OrcaOpt/orca5/orca
ulimit -s unlimited

#--------------------------------------
# Prepare Run
#--------------------------------------

export SLURMFILE=slurm-$SLURM_JOBID.out

#--------------------------------------
# Run CHARMM
#--------------------------------------

run orca file to replace
