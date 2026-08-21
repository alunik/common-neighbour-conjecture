#!/usr/bin/env python3
"""Multiply every Co1-orbit representative by eta in F_8.

The next two C++ stages identify the Co1 orbit containing each product.
"""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("triple_orbits", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with args.triple_orbits.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != TRIPLE_HEADER:
            raise RuntimeError("SCALAR QUERY FAILURE: bad triple header")
        rows = list(reader)
    if len(rows) != 9511:
        raise RuntimeError(
            f"SCALAR QUERY FAILURE: expected 9511 rows, got {len(rows)}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["query_id", "kind", "candidate_id", "a", "b", "c"])
        for orbit_index, row in enumerate(rows, 1):
            a, b, c = (int(row[name]) for name in ("a", "b", "c"))
            # In F_8 = F_2[eta]/(eta^3 + eta + 1),
            # eta * (a + eta*b + eta^2*c) = c + eta*(a+c) + eta^2*b.
            writer.writerow([orbit_index, "scalar", orbit_index, c, a ^ c, b])

    print(f"Scalar-action queries: {len(rows)}")


if __name__ == "__main__":
    main()
