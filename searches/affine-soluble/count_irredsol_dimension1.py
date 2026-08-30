#!/usr/bin/env python3
"""Count the dimension-one IRREDSOL actions in the census degree interval."""

from __future__ import annotations

import argparse
import json
from array import array


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lower", type=int, default=8192)
    parser.add_argument("--upper", type=int, default=2**24 - 1)
    args = parser.parse_args()
    if args.lower < 2 or args.upper < args.lower:
        raise SystemExit("invalid interval")

    # Smallest-prime-factor sieve.  For prime p, the dimension-one catalogue
    # has one cyclic subgroup for every divisor of p-1.
    spf = array("I", range(args.upper + 1))
    if args.upper >= 1:
        spf[1] = 1
    limit = int(args.upper**0.5)
    for p in range(2, limit + 1):
        if spf[p] != p:
            continue
        start = p * p
        for value in range(start, args.upper + 1, p):
            if spf[value] == value:
                spf[value] = p

    primes = 0
    actions = 0
    for p in range(max(2, args.lower), args.upper + 1):
        if spf[p] != p:
            continue
        primes += 1
        value = p - 1
        divisors = 1
        while value > 1:
            prime = spf[value]
            exponent = 0
            while value % prime == 0:
                value //= prime
                exponent += 1
            divisors *= exponent + 1
        actions += divisors

    print(
        json.dumps(
            {
                "schema": "IRREDSOL_DIMENSION1_CENSUS_V1",
                "lower_degree": args.lower,
                "upper_degree": args.upper,
                "prime_degrees": primes,
                "actions": actions,
                "base_size_1": primes,
                "base_size_2": actions - primes,
                "diameter_at_most_2": actions - primes,
                "diameter_at_least_3": 0,
                "proof": "nontrivial scalar subgroups act semiregularly on nonzero vectors",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
