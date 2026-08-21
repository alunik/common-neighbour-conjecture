#!/usr/bin/env python3
"""Check that every Co1 target orbit contains a sum of two regular vectors."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ORBIT_HEADER = [
    "orbit_index",
    "pair_index",
    "a",
    "b",
    "c",
    "first_orbit_size",
    "b_orbit_size",
    "c_orbit_size",
    "triple_stabilizer_order",
    "eta_orbit_index",
    "scalar_cycle_length",
    "maximal_regular",
]
CLASS_HEADER = [
    "query_id",
    "kind",
    "candidate_id",
    "pair_index",
    "c_rep",
    "triple_stabilizer_order",
]
CANDIDATE_HEADER = [
    "candidate_id",
    "target_orbit",
    "source_orbit",
    "scalar_power",
    "word_code",
    "x_a",
    "x_b",
    "x_c",
    "y_a",
    "y_b",
    "y_c",
    "z_a",
    "z_b",
    "z_c",
]
def fail(message: str) -> None:
    raise RuntimeError(f"WITNESS ANALYSIS FAILURE: {message}")


def read_orbits(path: Path) -> list[dict[str, int]]:
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ORBIT_HEADER:
            fail("bad classified-orbit header")
        rows = [{key: int(value) for key, value in row.items()} for row in reader]
    if len(rows) != 9511:
        fail(f"expected 9511 orbit rows, got {len(rows)}")
    if [row["orbit_index"] for row in rows] != list(range(1, 9512)):
        fail("noncontiguous orbit indices")
    return rows


def read_classifications(
    directory: Path, key_to_orbit: dict[tuple[int, int], int]
) -> dict[int, int]:
    result: dict[int, int] = {}
    for pair_index in range(1, 47):
        path = directory / f"pair-{pair_index}.tsv"
        with path.open(newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != CLASS_HEADER:
                fail(f"bad classification header for pair {pair_index}")
            for raw in reader:
                query_id = int(raw["query_id"])
                if raw["kind"] != "z":
                    fail(f"non-z query {query_id}")
                if int(raw["candidate_id"]) != query_id:
                    fail(f"candidate/query mismatch at {query_id}")
                if int(raw["pair_index"]) != pair_index:
                    fail(f"pair mismatch at query {query_id}")
                key = (pair_index, int(raw["c_rep"]))
                orbit_index = key_to_orbit.get(key)
                if orbit_index is None:
                    fail(f"unknown orbit key {key}")
                if (
                    int(raw["triple_stabilizer_order"])
                    != orbits_global[orbit_index - 1][
                        "triple_stabilizer_order"
                    ]
                ):
                    fail(f"stabilizer mismatch at query {query_id}")
                if query_id in result:
                    fail(f"duplicate query {query_id}")
                result[query_id] = orbit_index
    return result


orbits_global: list[dict[str, int]] = []


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("classified_orbits", type=Path)
    parser.add_argument("candidate_tsv", type=Path)
    parser.add_argument("classification_directory", type=Path)
    args = parser.parse_args()

    global orbits_global
    orbits_global = read_orbits(args.classified_orbits)
    key_to_orbit = {
        (row["pair_index"], row["c"]): row["orbit_index"]
        for row in orbits_global
    }
    if len(key_to_orbit) != len(orbits_global):
        fail("duplicate orbit key")
    classifications = read_classifications(
        args.classification_directory, key_to_orbit
    )

    selected: list[int | None] = [None] * (len(orbits_global) + 1)
    trial_counts = [0] * (len(orbits_global) + 1)
    candidate_count = 0
    with args.candidate_tsv.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != CANDIDATE_HEADER:
            fail("bad candidate header")
        for raw in reader:
            values = {key: int(value) for key, value in raw.items()}
            candidate_count += 1
            if values["candidate_id"] != candidate_count:
                fail("noncontiguous candidate id")
            target_index = values["target_orbit"]
            if not (1 <= target_index <= len(orbits_global)):
                fail(f"bad target at candidate {candidate_count}")
            trial_counts[target_index] += 1
            if selected[target_index] is not None:
                continue
            target = orbits_global[target_index - 1]
            source_index = values["source_orbit"]
            if not (1 <= source_index <= len(orbits_global)):
                fail(f"bad source at candidate {candidate_count}")
            source = orbits_global[source_index - 1]
            if not source["maximal_regular"]:
                fail(f"nonregular source at candidate {candidate_count}")
            if [values[f"x_{name}"] for name in ("a", "b", "c")] != [
                target[name] for name in ("a", "b", "c")
            ]:
                fail(f"target representative mismatch at {candidate_count}")
            for name in ("a", "b", "c"):
                if values[f"x_{name}"] ^ values[f"y_{name}"] != values[
                    f"z_{name}"
                ]:
                    fail(f"sum identity fails at candidate {candidate_count}")

            y_orbit = source_index
            for _ in range(values["scalar_power"]):
                y_orbit = orbits_global[y_orbit - 1]["eta_orbit_index"]
            if not orbits_global[y_orbit - 1]["maximal_regular"]:
                fail(f"generated y is not regular at {candidate_count}")
            z_orbit = classifications.get(candidate_count)
            if z_orbit is None:
                fail(f"missing classification at candidate {candidate_count}")
            if not orbits_global[z_orbit - 1]["maximal_regular"]:
                continue

            selected[target_index] = trial_counts[target_index]

    if set(classifications) != set(range(1, candidate_count + 1)):
        fail("classification query set is incomplete")
    missing = [index for index in range(1, len(selected)) if selected[index] is None]
    if missing:
        fail(
            f"{len(missing)} targets lack witnesses; first missing {missing[:20]}"
        )

    attempts = [int(selected[index]) for index in range(1, len(selected))]
    print(f"Target orbits checked: {len(orbits_global)}")
    print(
        f"Regular summands tried: {sum(attempts)} total, "
        f"at most {max(attempts)} for one target"
    )
    print("R + R = V")


if __name__ == "__main__":
    main()
