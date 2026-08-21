#!/usr/bin/env bash
#SBATCH --job-name=fi22-outer-data
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

gap_bin=${GAP_BIN:-gap}
"$gap_bin" -q -b sporadic/large/fi22/export_outer_actions.g
