#!/usr/bin/env python3
"""Recover the action of F_8^* on the Co1-orbits and find regular vectors."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


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
CLASS_HEADER = [
    "query_id",
    "kind",
    "candidate_id",
    "pair_index",
    "c_rep",
    "triple_stabilizer_order",
]


def fail(message: str) -> None:
    raise RuntimeError(f"SCALAR ANALYSIS FAILURE: {message}")


def read_rows(path: Path) -> list[dict[str, int]]:
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TRIPLE_HEADER:
            fail("bad triple header")
        rows = [{key: int(value) for key, value in row.items()} for row in reader]
    if len(rows) != 9511:
        fail(f"expected 9511 triple rows, got {len(rows)}")
    return rows


def read_classifications(directory: Path) -> dict[int, tuple[int, int, int]]:
    result: dict[int, tuple[int, int, int]] = {}
    for pair_index in range(1, 47):
        path = directory / f"pair-{pair_index}.tsv"
        with path.open(newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != CLASS_HEADER:
                fail(f"bad classification header for pair {pair_index}")
            for raw in reader:
                query_id = int(raw["query_id"])
                if raw["kind"] != "scalar":
                    fail(f"non-scalar query {query_id}")
                if int(raw["candidate_id"]) != query_id:
                    fail(f"candidate/query mismatch at {query_id}")
                if int(raw["pair_index"]) != pair_index:
                    fail(f"pair mismatch at query {query_id}")
                value = (
                    pair_index,
                    int(raw["c_rep"]),
                    int(raw["triple_stabilizer_order"]),
                )
                if query_id in result:
                    fail(f"duplicate query {query_id}")
                result[query_id] = value
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("triple_orbits", type=Path)
    parser.add_argument("classification_directory", type=Path)
    parser.add_argument("classified_output", type=Path)
    parser.add_argument("summary_output", type=Path)
    args = parser.parse_args()

    rows = read_rows(args.triple_orbits)
    classifications = read_classifications(args.classification_directory)
    if set(classifications) != set(range(1, len(rows) + 1)):
        fail("classification query set is incomplete")

    key_to_index: dict[tuple[int, int], int] = {}
    for orbit_index, row in enumerate(rows, 1):
        key = (row["pair_index"], row["c"])
        if key in key_to_index:
            fail(f"duplicate orbit key {key}")
        key_to_index[key] = orbit_index

    eta_map: list[int] = [0] * (len(rows) + 1)
    for orbit_index in range(1, len(rows) + 1):
        pair_index, c_rep, classified_stabilizer = classifications[orbit_index]
        target = key_to_index.get((pair_index, c_rep))
        if target is None:
            fail(f"unknown eta target at orbit {orbit_index}")
        if classified_stabilizer != rows[target - 1]["triple_stabilizer_order"]:
            fail(f"eta target stabilizer mismatch at orbit {orbit_index}")
        if (
            rows[orbit_index - 1]["triple_stabilizer_order"]
            != rows[target - 1]["triple_stabilizer_order"]
        ):
            fail(f"eta fails to preserve stabilizer at orbit {orbit_index}")
        eta_map[orbit_index] = target

    if set(eta_map[1:]) != set(range(1, len(rows) + 1)):
        fail("eta map is not a permutation")

    cycle_lengths = [0] * (len(rows) + 1)
    for orbit_index in range(1, len(rows) + 1):
        current = orbit_index
        length = 0
        while True:
            current = eta_map[current]
            length += 1
            if current == orbit_index:
                break
            if length > 7:
                fail(f"eta cycle longer than 7 at orbit {orbit_index}")
        if length not in (1, 7):
            fail(f"eta cycle has length {length} at orbit {orbit_index}")
        cycle_lengths[orbit_index] = length

    regular_co1 = sum(
        row["triple_stabilizer_order"] == 1 for row in rows
    )
    maximal_regular = [
        index
        for index, row in enumerate(rows, 1)
        if row["triple_stabilizer_order"] == 1 and cycle_lengths[index] == 7
    ]
    fixed_regular = regular_co1 - len(maximal_regular)
    if len(maximal_regular) % 7:
        fail("maximal regular Co1-orbit count is not divisible by 7")

    args.classified_output.parent.mkdir(parents=True, exist_ok=True)
    output_header = [
        "orbit_index",
        *TRIPLE_HEADER,
        "eta_orbit_index",
        "scalar_cycle_length",
        "maximal_regular",
    ]
    with args.classified_output.open(
        "w", newline="", encoding="ascii"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_header)
        for orbit_index, row in enumerate(rows, 1):
            writer.writerow(
                [
                    orbit_index,
                    *(row[name] for name in TRIPLE_HEADER),
                    eta_map[orbit_index],
                    cycle_lengths[orbit_index],
                    int(
                        row["triple_stabilizer_order"] == 1
                        and cycle_lengths[orbit_index] == 7
                    ),
                ]
            )

    fixed_cycles = sum(length == 1 for length in cycle_lengths[1:])
    seven_cycles = sum(length == 7 for length in cycle_lengths[1:]) // 7
    with args.summary_output.open(
        "w", newline="", encoding="ascii"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["quantity", "value"])
        writer.writerow(["co1_orbit_count", len(rows)])
        writer.writerow(["eta_fixed_co1_orbit_count", fixed_cycles])
        writer.writerow(["eta_seven_cycle_count", seven_cycles])
        writer.writerow(["regular_co1_orbit_count", regular_co1])
        writer.writerow(["scalar_fixed_regular_co1_orbit_count", fixed_regular])
        writer.writerow(["maximal_regular_co1_orbit_count", len(maximal_regular)])
        writer.writerow(
            ["maximal_regular_c7_co1_orbit_count", len(maximal_regular) // 7]
        )

    print(
        f"Scalar action: {fixed_cycles} fixed Co1-orbits and "
        f"{seven_cycles} cycles of length 7"
    )
    print(
        f"Regular orbits for F_8^* x Co1: {len(maximal_regular) // 7}"
    )


if __name__ == "__main__":
    main()
