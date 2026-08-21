#!/usr/bin/env bash
# Run the Co1 computation from the ATLAS generators to the final sumset test.
# All intermediate files live in a new temporary directory and are removed
# when the program exits.

set -euo pipefail

source_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
gap_bin=${GAP:-gap}
cxx=${CXX:-g++}
jobs=${JOBS:-4}
trials=${TRIALS:-128}

if ! command -v "$gap_bin" >/dev/null 2>&1; then
    echo "GAP was not found; set GAP to the GAP executable." >&2
    exit 1
fi
if ! command -v "$cxx" >/dev/null 2>&1; then
    echo "A C++ compiler was not found; set CXX to the compiler." >&2
    exit 1
fi

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/co1-f8.XXXXXX")
cleanup() {
    rm -rf -- "$work_dir"
}
trap cleanup EXIT

mkdir -p \
    "$work_dir/bin" \
    "$work_dir/pair_stabilizers" \
    "$work_dir/triple_orbits" \
    "$work_dir/scalars/normalized" \
    "$work_dir/scalars/classified" \
    "$work_dir/sums/normalized" \
    "$work_dir/sums/classified"

compile() {
    "$cxx" -O3 -std=c++17 -Wall -Wextra -Werror -pedantic \
        "$source_dir/$1.cpp" -o "$work_dir/bin/$1"
}

compile vector_orbits
compile pair_orbits
compile triple_orbits
compile normalize_queries
compile classify_queries
compile sum_queries

export CO1_WORK_DIR="$work_dir"
"$gap_bin" -q "$source_dir/group_data.g"
"$work_dir/bin/vector_orbits" \
    "$work_dir/group_data.txt" "$work_dir/vectors"
"$gap_bin" -q "$source_dir/pair_setup.g"
"$work_dir/bin/pair_orbits" \
    "$work_dir/pair_input.txt" "$work_dir/pair_orbits.tsv"

pair_count=$(($(wc -l < "$work_dir/pair_orbits.tsv") - 1))
if (( pair_count < 1 )); then
    echo "No pair orbits were found." >&2
    exit 1
fi

export source_dir work_dir gap_bin
seq 1 "$pair_count" | xargs -P "$jobs" -I '{}' bash -c '
    set -euo pipefail
    pair=$1
    CO1_PAIR_INDEX=$pair "$gap_bin" -q "$source_dir/pair_stabilizer.g"
    "$work_dir/bin/triple_orbits" \
        "$work_dir/pair_stabilizers/pair-$pair.txt" \
        "$work_dir/triple_orbits/pair-$pair.tsv"
' _ '{}'

python3 "$source_dir/merge_orbits.py" \
    "$work_dir/pair_orbits.tsv" \
    "$work_dir/triple_orbits" \
    "$work_dir/orbits.tsv" \
    "$work_dir/pair_summary.tsv"

# Determine the permutation of Co1-orbits induced by multiplication by eta.
python3 "$source_dir/scalar_queries.py" \
    "$work_dir/orbits.tsv" "$work_dir/scalars/queries.tsv"
"$work_dir/bin/normalize_queries" \
    "$work_dir/group_data.txt" \
    "$work_dir/pair_input.txt" \
    "$work_dir/pair_orbits.tsv" \
    "$work_dir/scalars/queries.tsv" \
    "$work_dir/scalars/normalized"

export phase=scalars
seq 1 "$pair_count" | xargs -P "$jobs" -I '{}' bash -c '
    set -euo pipefail
    pair=$1
    "$work_dir/bin/classify_queries" \
        "$work_dir/pair_stabilizers/pair-$pair.txt" \
        "$work_dir/triple_orbits/pair-$pair.tsv" \
        "$work_dir/$phase/normalized/pair-$pair.tsv" \
        "$work_dir/$phase/classified/pair-$pair.tsv"
' _ '{}'

python3 "$source_dir/scalar_orbits.py" \
    "$work_dir/orbits.tsv" \
    "$work_dir/scalars/classified" \
    "$work_dir/regular_orbits.tsv" \
    "$work_dir/scalar_summary.tsv"

# Try regular summands, classify their complements, and check every target.
"$work_dir/bin/sum_queries" \
    "$work_dir/group_data.txt" \
    "$work_dir/regular_orbits.tsv" \
    "$trials" \
    "$work_dir/sums/queries.tsv" \
    "$work_dir/sums/candidates.tsv"
"$work_dir/bin/normalize_queries" \
    "$work_dir/group_data.txt" \
    "$work_dir/pair_input.txt" \
    "$work_dir/pair_orbits.tsv" \
    "$work_dir/sums/queries.tsv" \
    "$work_dir/sums/normalized"

export phase=sums
seq 1 "$pair_count" | xargs -P "$jobs" -I '{}' bash -c '
    set -euo pipefail
    pair=$1
    "$work_dir/bin/classify_queries" \
        "$work_dir/pair_stabilizers/pair-$pair.txt" \
        "$work_dir/triple_orbits/pair-$pair.tsv" \
        "$work_dir/$phase/normalized/pair-$pair.tsv" \
        "$work_dir/$phase/classified/pair-$pair.tsv"
' _ '{}'

python3 "$source_dir/check_sums.py" \
    "$work_dir/regular_orbits.tsv" \
    "$work_dir/sums/candidates.tsv" \
    "$work_dir/sums/classified"
