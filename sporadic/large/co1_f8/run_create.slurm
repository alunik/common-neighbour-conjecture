#!/usr/bin/env bash
#SBATCH --job-name=co1-f8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G

set -euo pipefail

case_dir=${CASE_DIR:-${SLURM_SUBMIT_DIR:-}}
if [[ ! -x "$case_dir/run.sh" ]]; then
    echo "Submit this file from its own directory, or set CASE_DIR." >&2
    exit 1
fi
export JOBS=${SLURM_CPUS_PER_TASK:-12}
"$case_dir/run.sh"
