#!/usr/bin/env bash
#SBATCH --job-name=fi22-summary
#SBATCH --partition=cpu
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

out=sporadic/large/fi22/generated
parts=$out/checked_parts
result=$out/checked_sums.tsv

{
  head -n 1 "$parts/checked_0.tsv"
  for part in $(seq 0 95); do
    tail -n +2 "$parts/checked_${part}.tsv"
  done
} > "$result"

awk -F '\t' '
  NR == 1 { next }
  {
    expected = NR - 2;
    if ($1 != expected || $10 != 1) exit 1;
  }
  END {
    if (NR != 6107) exit 1;
  }
' "$result"

cat <<'EOF'
5,363 inner targets and 743 outer targets were considered.
Each of the 6,106 nonregular target representatives is a sum of two
vectors with trivial stabiliser in 3.Fi22.2.

Together with the immediate decompositions for zero and for regular
targets, this gives R_H + R_H = V for H = 3.Fi22.2 on GF(2)^54.
The same equality holds for 3.Fi22 on GF(2)^54 and, through restriction
of scalars, for 3.Fi22 on GF(4)^27.
EOF
