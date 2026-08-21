#!/usr/bin/env bash
set -euo pipefail

# Co3 on GF(4)^22 is represented as ordered pairs of binary vectors.
# GAP exports the standard generators and five first-coordinate stabilisers;
# the C++ programs then enumerate the binary and ordered-pair orbits.

here="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/co3-f4.XXXXXX")"
trap 'rm -rf -- "$work"' EXIT

gap_bin="${GAP_BIN:-gap}"
cxx_bin="${CXX:-c++}"

"$gap_bin" -q -b "$here/module.g" | tr -d '\r' > "$work/module.txt"
"$gap_bin" -q -b "$here/stabilisers.g" | tr -d '\r' > "$work/stabilisers.txt"

"$cxx_bin" -O3 -std=c++17 -Wall -Wextra -Wpedantic \
  "$here/orbits.cpp" -o "$work/orbits"
"$cxx_bin" -O3 -std=c++17 -Wall -Wextra -Wpedantic \
  "$here/pairs.cpp" -o "$work/pairs"
"$cxx_bin" -O3 -std=c++17 -Wall -Wextra -Wpedantic \
  "$here/sumset.cpp" -o "$work/sumset"

"$work/orbits" "$work/module.txt" > "$work/orbits.txt"
"$work/pairs" "$work/stabilisers.txt" > "$work/pairs.txt"
"$work/sumset" "$work/module.txt" "$work/stabilisers.txt" \
  > "$work/sumset.txt"

grep -E '^(orbit_count|covered|total_pair_orbits|total_regular_pair_orbits|regular_Co3_vector_orbits|regular_Co3_vectors|regular_full_scalar_vectors|target_orbits|all_target_orbits_have_verified_decomposition|sumset_equals_GF4_22)' \
  "$work/orbits.txt" "$work/pairs.txt" "$work/sumset.txt"
