#!/bin/bash
#SBATCH -p multicore
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --time=10:00:00

# Your commands here
bash run_all_experiments.sh
