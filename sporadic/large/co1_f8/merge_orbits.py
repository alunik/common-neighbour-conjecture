#!/usr/bin/env python3
"""Merge the 46 third-coordinate searches into the Co1 orbit list on W^3.

The script also checks every orbit-stabilizer identity and that the resulting
orbits have total size 2^72.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


CO1_ORDER = 4_157_776_806_543_360_000
W_SIZE = 1 << 24
PAIR_HEADER = [
    "pair_index",
    "a",
    "b",
    "first_orbit_size",
    "b_orbit_size",
    "pair_stabilizer_order",
]
TRIPLE_HEADER = [
    "pair_index",
    "a",
    "b",
    "c",
    "first_orbit_size",
    "b_orbit_size",
    "c_orbit_size",
    "triple_stabilizer_order",
]


def fail(message: str) -> None:
    raise RuntimeError(f"HIERARCHY AGGREGATION FAILURE: {message}")


def read_pair_rows(path: Path) -> list[dict[str, int]]:
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != PAIR_HEADER:
            fail("pair-orbit header")
        rows = [{key: int(value) for key, value in row.items()} for row in reader]
    if len(rows) != 46:
        fail(f"expected 46 pair rows, got {len(rows)}")
    if [row["pair_index"] for row in rows] != list(range(1, 47)):
        fail("pair indices are not contiguous")
    pair_mass = 0
    for row in rows:
        if (
            row["first_orbit_size"]
            * row["b_orbit_size"]
            * row["pair_stabilizer_order"]
            != CO1_ORDER
        ):
            fail(f"pair orbit-stabilizer identity at {row['pair_index']}")
        pair_mass += row["first_orbit_size"] * row["b_orbit_size"]
    if pair_mass != 1 << 48:
        fail(f"pair mass {pair_mass} != 2^48")
    return rows


def aggregate(
    pair_path: Path,
    triple_dir: Path,
    merged_path: Path,
    summary_path: Path,
) -> None:
    pair_rows = read_pair_rows(pair_path)
    global_mass = 0
    global_orbits = 0
    global_regular_orbits = 0

    with (
        merged_path.open("w", newline="", encoding="ascii") as merged_handle,
        summary_path.open("w", newline="", encoding="ascii") as summary_handle,
    ):
        merged = csv.writer(merged_handle, delimiter="\t", lineterminator="\n")
        summary = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        merged.writerow(TRIPLE_HEADER)
        summary.writerow(
            [
                "pair_index",
                "a",
                "b",
                "pair_orbit_size",
                "pair_stabilizer_order",
                "triple_orbit_count",
                "regular_triple_orbit_count",
                "c_mass",
            ]
        )

        for pair in pair_rows:
            pair_index = pair["pair_index"]
            triple_path = triple_dir / f"pair-{pair_index}.tsv"
            if not triple_path.is_file():
                fail(f"missing {triple_path}")
            c_mass = 0
            orbit_count = 0
            regular_count = 0
            previous_c = -1
            with triple_path.open(newline="", encoding="ascii") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != TRIPLE_HEADER:
                    fail(f"triple header at pair {pair_index}")
                for raw in reader:
                    row = {key: int(value) for key, value in raw.items()}
                    if (
                        row["pair_index"] != pair_index
                        or row["a"] != pair["a"]
                        or row["b"] != pair["b"]
                        or row["first_orbit_size"] != pair["first_orbit_size"]
                        or row["b_orbit_size"] != pair["b_orbit_size"]
                    ):
                        fail(f"triple/pair alignment at pair {pair_index}")
                    if not (0 <= row["c"] < W_SIZE) or row["c"] <= previous_c:
                        fail(f"non-increasing c representative at pair {pair_index}")
                    previous_c = row["c"]
                    c_orbit = row["c_orbit_size"]
                    triple_stabilizer = row["triple_stabilizer_order"]
                    if c_orbit <= 0 or triple_stabilizer <= 0:
                        fail(f"nonpositive orbit data at pair {pair_index}")
                    if c_orbit * triple_stabilizer != pair["pair_stabilizer_order"]:
                        fail(f"triple orbit-stabilizer identity at pair {pair_index}")
                    triple_mass = (
                        pair["first_orbit_size"]
                        * pair["b_orbit_size"]
                        * c_orbit
                    )
                    if triple_mass * triple_stabilizer != CO1_ORDER:
                        fail(f"global orbit-stabilizer identity at pair {pair_index}")
                    c_mass += c_orbit
                    global_mass += triple_mass
                    orbit_count += 1
                    global_orbits += 1
                    if triple_stabilizer == 1:
                        regular_count += 1
                        global_regular_orbits += 1
                    merged.writerow(row[key] for key in TRIPLE_HEADER)
            if c_mass != W_SIZE:
                fail(f"c mass at pair {pair_index}: {c_mass} != 2^24")
            summary.writerow(
                [
                    pair_index,
                    pair["a"],
                    pair["b"],
                    pair["first_orbit_size"] * pair["b_orbit_size"],
                    pair["pair_stabilizer_order"],
                    orbit_count,
                    regular_count,
                    c_mass,
                ]
            )

    if global_mass != 1 << 72:
        fail(f"global mass {global_mass} != 2^72")
    print(f"Co1 pair orbits on (F_2^24)^2: {len(pair_rows)}")
    print(
        "Co1 orbits on (F_2^24)^3: "
        f"{global_orbits} (regular {global_regular_orbits})"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pair_orbits", type=Path)
    parser.add_argument("triple_directory", type=Path)
    parser.add_argument("merged_output", type=Path)
    parser.add_argument("summary_output", type=Path)
    args = parser.parse_args()
    aggregate(
        args.pair_orbits,
        args.triple_directory,
        args.merged_output,
        args.summary_output,
    )


if __name__ == "__main__":
    main()
