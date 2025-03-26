#!/bin/bash
#SBATCH --job-name=bench
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu06

# Load modules
module load charmm/gcc-14.2-ompi-5.0-cuda-12.6-omm

# Make scratch working directory
cdir=`pwd`
wdir=/scratch/$USER.$SLURM_JOB_ID
if [ -d "$wdir" ]; then
   echo "Working directory already exists"
else
   echo "Make working directory $wdir"
   mkdir $wdir
fi

# Copy files to working directory
cp -rf $cdir/* $wdir

# Switch to working directory
cd $wdir

# Run job
charmm -i gpu_step5_production.inp -o gpu_step5_production.out
echo "Job done"

# Copy files back to current directory
cp -rf $wdir/* $cdir/.
echo "Files copied back to $cdir"
