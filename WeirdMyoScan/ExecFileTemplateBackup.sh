#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=1000

#--------------------------------------
# Modules
#--------------------------------------

module load intel/2019-compiler-intel-openmpi-4.1.2
my_orca=/home/kaeserj/data/GaussianOpt/OrcaOpt/orca5
ulimit -s 10420

#--------------------------------------
# Prepare Run
#--------------------------------------

export SLURMFILE=slurm-$SLURM_JOBID.out

#--------------------------------------
# Run CHARMM
#--------------------------------------

srun $my_orca input_file_path > output_file_path
