#!/usr/bin/env bash
#SBATCH --job-name=fi22-find-sums
#SBATCH --partition=cpu
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=24G
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

src=sporadic/large/fi22
out=$src/generated
threads=${SLURM_CPUS_PER_TASK:-64}

# Targets fixed by an element in the index-two subgroup.
build/fi22/bin/find_inner_sums \
  "$src" "$out/inner_sum_candidates.tsv" "$threads" 12000 40

awk '
  !/^#/ { count[$1]++ }
  END {
    for (i = 0; i < 5363; ++i)
      if (count[i] != 40) exit 1;
    print "Found 40 candidate sums for each of 5,363 inner targets.";
  }
' "$out/inner_sum_candidates.tsv"

# Targets fixed by an outer involution.  First reconcile the independently
# constructed GF(4)-to-GF(2) coordinate maps, then search for sums.
build/fi22/bin/combine_targets \
  "$out/outer_fixed/orbits_outer_1.txt" \
  "$out/outer_fixed/orbits_outer_2.txt" \
  "$out/outer_fixed/orbits_outer_3.txt" \
  "$out/inner_sum_candidates.tsv" \
  "$out/outer_fixed/outer_targets_unique.tsv"

build/fi22/bin/map_outer_targets \
  "$out/outer_fixed/outer_targets_unique.tsv" \
  "$out/fischer/f4_to_h54.tsv" \
  "$out/inner_sum_candidates.tsv" \
  "$out/outer_fixed/outer_targets_wl.tsv"

# The outer search keeps more alternatives because one target needs a later
# random-walk image than the others.
build/fi22/bin/find_outer_sums \
  "$src" "$out/outer_sum_candidates.tsv" "$threads" 60000 200

awk '
  !/^#/ { count[$1]++ }
  END {
    for (i = 5363; i < 6106; ++i)
      if (count[i] != 200) exit 1;
    print "Found 200 candidate sums for each of 743 outer targets.";
  }
' "$out/outer_sum_candidates.tsv"

{
  echo '# Fi22 candidate sums'
  echo '# tid class block coordinates target choice epoch attempt global_attempt norm zeros word first second'
  grep -v '^#' "$out/inner_sum_candidates.tsv"
  grep -v '^#' "$out/outer_sum_candidates.tsv"
} > "$out/all_sum_candidates.tsv"

FI22_FISCHER_DATA="$out/fischer" \
  build/fi22/bin/refine_candidates \
  "$out/all_sum_candidates.tsv" \
  "$out/refined_candidates.tsv"

echo "Refined the candidates for all 6,106 target representatives."
