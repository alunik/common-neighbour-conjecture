#!/usr/bin/env bash
#SBATCH --job-name=fi22-outer-orbits
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --array=0-2
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

src=sporadic/large/fi22
class=$((SLURM_ARRAY_TASK_ID + 1))

build/fi22/bin/enumerate_outer_orbits \
  "$src/generated/outer_fixed/outer_${class}.txt" \
  "$src/generated/outer_fixed/orbits_outer_${class}.txt" \
  "${SLURM_CPUS_PER_TASK:-32}"

tail -n 1 "$src/generated/outer_fixed/orbits_outer_${class}.txt"
