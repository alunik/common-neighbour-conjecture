#!/usr/bin/env bash
set -euo pipefail

# This is an exhaustive run (roughly ten minutes on one CREATE CPU core).
# All intermediate output lives in a temporary directory.

here="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
work="$(mktemp -d "${TMPDIR:-/tmp}/co1-f9.XXXXXX")"
trap 'rm -rf -- "$work"' EXIT

cxx_bin="${CXX:-c++}"
gap_bin="${GAP_BIN:-gap}"
pool_size="${POOL_SIZE:-4000}"

"$cxx_bin" -O3 -std=c++17 -Wall -Wextra -Wpedantic \
  "$here/search.cpp" -o "$work/search"

"$work/search" "$here/representation.txt" "$pool_size" \
  > "$work/candidates.tsv"

export CO1_F9_DIR="$here"
export CO1_F9_CANDIDATES="$work/candidates.tsv"
"$gap_bin" -q -b "$here/check.g"
