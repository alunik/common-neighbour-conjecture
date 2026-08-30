#!/usr/bin/env python3
"""List prime-field affine degrees not covered by the degree-8191 catalogue."""

from __future__ import annotations

import argparse
import math


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, math.isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [value for value, flag in enumerate(sieve) if flag]


def target_pairs(min_degree: int, max_degree: int, max_dimension: int) -> list[tuple[int, int, int]]:
    rows: list[tuple[int, int, int]] = []
    for dimension in range(2, max_dimension + 1):
        prime_limit = int(max_degree ** (1 / dimension)) + 2
        for prime in primes_up_to(prime_limit):
            degree = prime**dimension
            if min_degree <= degree <= max_degree:
                rows.append((dimension, prime, degree))
    return sorted(rows, key=lambda row: (row[0], row[1]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-degree", type=int, default=8192)
    parser.add_argument("--max-degree", type=int, default=1_000_000)
    parser.add_argument("--max-dimension", type=int, default=19)
    args = parser.parse_args()
    if args.min_degree < 2 or args.min_degree > args.max_degree:
        raise SystemExit("invalid degree interval")
    if args.max_dimension < 2:
        raise SystemExit("max dimension must be at least 2")

    rows = target_pairs(args.min_degree, args.max_degree, args.max_dimension)
    print("dimension\tprime\tdegree\tmethod")
    for dimension, prime, degree in rows:
        method = "IrreducibleSubgroups" if dimension in (2, 3) else "ClassicalMaximals_frontier"
        print(f"{dimension}\t{prime}\t{degree}\t{method}")

    direct = sum(dimension in (2, 3) for dimension, _, _ in rows)
    frontier = len(rows) - direct
    print(
        f"TARGET_SUMMARY total={len(rows)} direct={direct} frontier={frontier} "
        f"min_degree={args.min_degree} max_degree={args.max_degree} "
        f"max_dimension={args.max_dimension}",
        file=__import__("sys").stderr,
    )


if __name__ == "__main__":
    main()
