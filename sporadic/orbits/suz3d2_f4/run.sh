#!/usr/bin/env bash
set -euo pipefail

# The 24-dimensional GF(4)-module is represented as two binary components.
# GAP supplies the standard matrices and first-coordinate stabilisers; the
# packed C++ routines perform the exhaustive orbit and sumset calculations.

here="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/suz3d2-f4.XXXXXX")"
trap 'rm -rf -- "$work"' EXIT

gap_bin="${GAP_BIN:-gap}"
cxx_bin="${CXX:-c++}"

"$gap_bin" -q -b "$here/module.g" | tr -d '\r' > "$work/module.txt"
"$gap_bin" -q -b "$here/stabilisers.g" | tr -d '\r' > "$work/stabilisers.txt"

"$cxx_bin" -O3 -std=c++20 -Wall -Wextra -Wpedantic \
  "$here/orbits.cpp" -o "$work/orbits"
"$cxx_bin" -O3 -std=c++20 -Wall -Wextra -Wpedantic \
  "$here/pairs.cpp" -o "$work/pairs"
"$cxx_bin" -O3 -std=c++20 -Wall -Wextra -Wpedantic \
  "$here/sumset.cpp" -o "$work/sumset"

"$work/orbits" "$work/module.txt" > "$work/orbits.txt"
"$work/pairs" "$work/stabilisers.txt" "$work/pair-orbits.tsv" \
  > "$work/pairs.txt"
"$work/sumset" "$work/module.txt" "$work/stabilisers.txt" \
  "$work/sums.tsv" > "$work/sumset.txt"

grep -E '^(orbit_count|covered|total_pair_orbits|total_regular_pair_orbits|GF4_vector_orbits|regular_G_vector_orbits|regular_G_vectors|regular_Hmax_vector_orbits|regular_Hmax_vectors|G_sumset_equals_GF4_24|Hmax_sumset_equals_GF4_24)' \
  "$work/orbits.txt" "$work/pairs.txt" "$work/sumset.txt"
