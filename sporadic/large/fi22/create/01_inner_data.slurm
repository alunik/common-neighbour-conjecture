#!/usr/bin/env bash
#SBATCH --job-name=fi22-inner-data
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

src=sporadic/large/fi22
mkdir -p "$src/generated"
gap_bin=${GAP_BIN:-gap}

"$gap_bin" -q -b "$src/export_core_geometry.g"
"$gap_bin" -q -b "$src/export_involution_spaces.g"
"$gap_bin" -q -b "$src/export_prime_order_spaces.g"

echo "Wrote the prime-order fixed-space actions."
