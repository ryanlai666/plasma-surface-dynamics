#!/bin/bash
#SBATCH --job-name=plasma-surface
#SBATCH --array=0-15
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%A_%a.out
set -euo pipefail
# Submit from the repository root after installing into your cluster environment.
# Site-specific account, partition, and module/environment setup must be added.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
python -m plasma_surface.cli sweep --samples 4096 --seed 42 \
  --shard "$SLURM_ARRAY_TASK_ID" --shards 16 \
  --output "outputs/hpc/shard_${SLURM_ARRAY_TASK_ID}.csv"
