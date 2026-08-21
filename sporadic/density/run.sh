#!/usr/bin/env bash
set -euo pipefail

here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build="$(mktemp -d)"
trap 'rm -rf -- "${build}"' EXIT

compiler="${CXX:-c++}"

"${compiler}" -std=c++17 -O2 -Wall -Wextra -pedantic \
  "${here}/co3_f5.cpp" -o "${build}/co3_f5"
"${compiler}" -std=c++17 -O2 -Wall -Wextra -pedantic \
  "${here}/suz3_f16.cpp" -o "${build}/suz3_f16"

"${build}/co3_f5"
printf '\n'
"${build}/suz3_f16"

if [[ "${1:-}" == "--full" ]]; then
  gap_bin="${GAP_BIN:-gap}"
  printf '\nRunning the GAP computations...\n\n'
  "${gap_bin}" -q -b "${here}/co3_f5.g"
  printf '\n'
  "${gap_bin}" -q -b "${here}/suz3_f16.g"
elif [[ $# -ne 0 ]]; then
  printf 'usage: %s [--full]\n' "$0" >&2
  exit 2
fi
