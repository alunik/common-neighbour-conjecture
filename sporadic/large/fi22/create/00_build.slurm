#!/usr/bin/env bash
#SBATCH --job-name=fi22-build
#SBATCH --partition=cpu
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=build/fi22/logs/%x-%j.out
#SBATCH --error=build/fi22/logs/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"

src=sporadic/large/fi22
bin=build/fi22/bin
mkdir -p "$bin"

g++ -std=c++20 -O3 -pthread -Wall -Wextra -Wpedantic \
  "$src/enumerate_inner_orbits.cpp" -o "$bin/enumerate_inner_orbits"
g++ -std=c++20 -O3 -pthread -Wall -Wextra -Wpedantic \
  "$src/enumerate_outer_orbits.cpp" -o "$bin/enumerate_outer_orbits"
g++ -std=c++20 -O3 -pthread -Wall -Wextra -Wpedantic \
  "$src/find_inner_sums.cpp" -o "$bin/find_inner_sums"
g++ -std=c++20 -O3 -pthread -Wall -Wextra -Wpedantic \
  "$src/find_outer_sums.cpp" -o "$bin/find_outer_sums"
g++ -std=c++20 -O3 -Wall -Wextra -Wpedantic \
  "$src/combine_targets.cpp" -o "$bin/combine_targets"
g++ -std=c++20 -O3 -Wall -Wextra -Wpedantic \
  "$src/map_outer_targets.cpp" -o "$bin/map_outer_targets"
g++ -std=c++20 -O3 -pthread -Wall -Wextra -Wpedantic \
  "$src/refine_candidates.cpp" -o "$bin/refine_candidates"

echo "Built the Fi22 programs."
