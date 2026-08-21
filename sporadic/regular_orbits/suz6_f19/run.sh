#!/usr/bin/env bash
set -euo pipefail

source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/suz6-f19.XXXXXX")"
trap 'rm -rf -- "$work_dir"' EXIT

cxx="${CXX:-g++}"
gap="${GAP_BIN:-gap}"

"$cxx" -O3 -std=c++17 -Wall -Wextra -pedantic \
  "$source_dir/main.cpp" -o "$work_dir/orbit_search"

"$work_dir/orbit_search" "$work_dir/A" A
"$work_dir/orbit_search" "$work_dir/B" B

export SUZ6_F19_SOURCE="$source_dir"
export SUZ6_F19_WORK="$work_dir"
"$gap" -q "$source_dir/compute.g"
