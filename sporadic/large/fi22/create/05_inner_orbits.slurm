#!/usr/bin/env bash
#SBATCH --job-name=fi22-inner-orbits
#SBATCH --partition=cpu
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --array=0-7%3
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

src=sporadic/large/fi22
classes=(2A 2B 2C 3A 3B 3C 3D 5A)
class=${classes[${SLURM_ARRAY_TASK_ID:?}]}

build/fi22/bin/enumerate_inner_orbits \
  "$src/generated/cls_${class}.txt" \
  "$src/generated/orbits_${class}.txt" \
  "${SLURM_CPUS_PER_TASK:-32}"

tail -n 1 "$src/generated/orbits_${class}.txt"
