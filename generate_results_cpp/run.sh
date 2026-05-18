#!/bin/bash
#SBATCH --job-name=scite_rna
#SBATCH --output=logs/output_%A_%a.log
#SBATCH --error=logs/error_%A_%a.err
#SBATCH --array=0-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00

INPUT_DIR="50c500m"

cd build

./SCITE-RNA "$INPUT_DIR" "$SLURM_ARRAY_TASK_ID"
