#!/usr/bin/env bash
#SBATCH --job-name=fi22-check-sums
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=0-95%32
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

src=sporadic/large/fi22
out=$src/generated
part=$out/checked_parts
mkdir -p "$part"

chunk=64
start=$((SLURM_ARRAY_TASK_ID * chunk))
end=$((start + chunk - 1))
(( end > 6105 )) && end=6105

input="$part/input_${SLURM_ARRAY_TASK_ID}.tsv"
output="$part/checked_${SLURM_ARRAY_TASK_ID}.tsv"
{
  echo '# Fi22 refined candidates for one target range'
  awk -v first="$start" -v last="$end" \
    '!/^#/ && $1 >= first && $1 <= last { print }' \
    "$out/refined_candidates.tsv"
} > "$input"

export FI22_REFINED_INPUT="$input"
export FI22_STABILISER_OUTPUT="$output"
export FI22_CONTAINER_LIMIT=10000
gap_bin=${GAP_BIN:-gap}
"$gap_bin" -q -b "$src/check_second_summands.g"

expected=$((end - start + 1))
actual=$(($(wc -l < "$output") - 1))
[[ "$actual" -eq "$expected" ]]
echo "Checked targets $start through $end."
