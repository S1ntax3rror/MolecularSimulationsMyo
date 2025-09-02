#!/bin/bash
#SBATCH --job-name=Analysis_Main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --array=0-50

echo $HOSTNAME

mapfile -t arr < subfolders.txt

# Access by index
echo "${arr[$SLURM_ARRAY_TASK_ID]}"   # first line

# Cast whole array to string
str="${arr[$SLURM_ARRAY_TASK_ID]}"
echo "$str"
cd $str

python3 plot_plus_gen_npz.py

echo "All done!"
