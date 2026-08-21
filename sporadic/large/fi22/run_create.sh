#!/usr/bin/env bash

# Submit the Fi22 calculation from the root of the repository.  Slurm
# dependencies keep each generated input ahead of the program that reads it.

set -euo pipefail

repo=$(cd "$(dirname "$0")/../../.." && pwd)
cd "$repo"
mkdir -p build/fi22/logs

create=sporadic/large/fi22/create
build_job=$(sbatch --parsable "$create/00_build.slurm")

inner_data=$(sbatch --parsable --dependency="afterok:$build_job" \
  "$create/01_inner_data.slurm")
regular_seed=$(sbatch --parsable --dependency="afterok:$inner_data" \
  "$create/02_regular_seed.slurm")
outer_data=$(sbatch --parsable --dependency="afterok:$build_job" \
  "$create/03_outer_data.slurm")
fischer_data=$(sbatch --parsable --dependency="afterok:$build_job" \
  "$create/04_fischer_data.slurm")

inner_orbits=$(sbatch --parsable \
  --dependency="afterok:$inner_data" "$create/05_inner_orbits.slurm")
outer_orbits=$(sbatch --parsable \
  --dependency="afterok:$outer_data" "$create/06_outer_orbits.slurm")

find_sums=$(sbatch --parsable \
  --dependency="afterok:$regular_seed:$fischer_data:$inner_orbits:$outer_orbits" \
  "$create/07_find_sums.slurm")
check_sums=$(sbatch --parsable \
  --dependency="afterok:$find_sums" "$create/08_check_sums.slurm")
summary=$(sbatch --parsable \
  --dependency="afterok:$check_sums" "$create/09_summary.slurm")

echo "Submitted Fi22 calculation.  Final job: $summary"
