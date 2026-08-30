#!/usr/bin/env python3
"""Exact distinguishing-colouring reductions for every CD18 wreath overgroup.

For a component group L on Gamma with point stabiliser H and r regular
H-orbits, the neighbours of a fixed point in L wr P are counted by

    d_P(r) * |H|^m,

where d_P(r) is the number of distinguishing r-colourings of the m points.
Here m <= 5, so d_P is computed exactly from all set partitions.  A strict
half-density bound closes the full wreath overgroup and every subgroup in it.

There is also a top-group-independent route.  If r regular H-orbits are
guaranteed, then for every target point at least

    ceil((2*r*|H| - |Gamma|) / |H|)

regular colours at the base point also occur regularly at the target.  If
this is at least m, Hall's condition gives m distinct such colours for the m
coordinates.  The resulting all-singleton colouring distinguishes every
permutation group of degree m, so it gives a common neighbour for every
subgroup of the full wreath product.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable


class AnalysisError(RuntimeError):
    pass


def require(value: bool, message: str) -> None:
    if not value:
        raise AnalysisError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def create(path: Path, content: bytes, mode: int = 0o400) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, mode)
    try:
        os.write(descriptor, content)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fresh(path: Path) -> Path:
    path = path.resolve()
    require(path != Path("/") and path != Path.home() and
            not path.exists() and not path.is_symlink(), "output is not fresh/safe")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.mkdir(mode=0o700)
    return path


Permutation = tuple[int, ...]


def perm(*cycles: tuple[int, ...], degree: int) -> Permutation:
    result = list(range(degree))
    for cycle in cycles:
        for source, target in zip(cycle, cycle[1:] + cycle[:1]):
            result[source] = target
    return tuple(result)


def compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(right[left[index]] for index in range(len(left)))


def closure(degree: int, generators: Iterable[Permutation]) -> tuple[Permutation, ...]:
    identity = tuple(range(degree))
    group = {identity}
    queue = [identity]
    generators = tuple(generators)
    while queue:
        current = queue.pop()
        for generator in generators:
            product = compose(current, generator)
            if product not in group:
                group.add(product)
                queue.append(product)
    return tuple(sorted(group))


def top_groups() -> dict[int, tuple[tuple[str, tuple[Permutation, ...], int], ...]]:
    raw: dict[int, tuple[tuple[str, tuple[Permutation, ...], int], ...]] = {
        2: (("S2", (perm((0, 1), degree=2),), 2),),
        3: (
            ("C3", (perm((0, 1, 2), degree=3),), 3),
            ("S3", (perm((0, 1, 2), degree=3), perm((0, 1), degree=3)), 6),
        ),
        4: (
            ("C4", (perm((0, 1, 2, 3), degree=4),), 4),
            ("V4", (perm((0, 1), (2, 3), degree=4),
                    perm((0, 2), (1, 3), degree=4)), 4),
            ("D8", (perm((0, 1, 2, 3), degree=4), perm((1, 3), degree=4)), 8),
            ("A4", (perm((0, 1, 2), degree=4), perm((0, 1, 3), degree=4)), 12),
            ("S4", (perm((0, 1, 2, 3), degree=4), perm((0, 1), degree=4)), 24),
        ),
        5: (
            ("C5", (perm((0, 1, 2, 3, 4), degree=5),), 5),
            ("D10", (perm((0, 1, 2, 3, 4), degree=5),
                     perm((1, 4), (2, 3), degree=5)), 10),
            ("F20", (perm((0, 1, 2, 3, 4), degree=5),
                     perm((1, 2, 4, 3), degree=5)), 20),
            ("A5", (perm((0, 1, 2, 3, 4), degree=5),
                    perm((0, 1, 2), degree=5)), 60),
            ("S5", (perm((0, 1, 2, 3, 4), degree=5),
                    perm((0, 1), degree=5)), 120),
        ),
    }
    answer: dict[int, tuple[tuple[str, tuple[Permutation, ...], int], ...]] = {}
    for degree, rows in raw.items():
        checked = []
        for name, generators, expected_order in rows:
            group = closure(degree, generators)
            require(len(group) == expected_order, f"bad top group order: {name}")
            require({g[0] for g in group} == set(range(degree)),
                    f"intransitive top group: {name}")
            checked.append((name, group, expected_order))
        answer[degree] = tuple(checked)
    require(sum(len(rows) for rows in answer.values()) == 13, "top group census mismatch")
    return answer


def set_partitions(degree: int) -> Iterable[tuple[int, ...]]:
    # Restricted-growth strings, hence one representative per set partition.
    def extend(prefix: tuple[int, ...]) -> Iterable[tuple[int, ...]]:
        if len(prefix) == degree:
            yield prefix
            return
        for value in range(max(prefix) + 2):
            yield from extend(prefix + (value,))

    yield from extend((0,))


def stabilises_partition(action: Permutation, partition: tuple[int, ...]) -> bool:
    return all(partition[action[index]] == partition[index]
               for index in range(len(partition)))


def distinguishing_coefficients(group: tuple[Permutation, ...]) -> dict[int, int]:
    identity = tuple(range(len(group[0])))
    coefficients: dict[int, int] = {}
    for partition in set_partitions(len(identity)):
        stabiliser = {g for g in group if stabilises_partition(g, partition)}
        if stabiliser == {identity}:
            blocks = max(partition) + 1
            coefficients[blocks] = coefficients.get(blocks, 0) + 1
    return coefficients


def falling(value: int, length: int) -> int:
    if value < length:
        return 0
    return math.prod(range(value - length + 1, value + 1))


def distinguishing_colourings(coefficients: dict[int, int], colours: int) -> int:
    return sum(count * falling(colours, blocks)
               for blocks, count in coefficients.items())


def ceil_fraction(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--component-results", type=Path, required=True)
    parser.add_argument("--component-audit", type=Path, required=True)
    parser.add_argument("--arithmetic-inventory", type=Path, required=True)
    parser.add_argument("--component-inventory", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        audit = json.loads(args.component_audit.read_text())
        require(audit.get("schema") in {"CD18_COMPONENT_FPR_AUDIT_V1",
                                        "CD18_COMPONENT_PAIR_FPR_AUDIT_V1"} and
                audit.get("complete") is True and audit.get("task_count") == 381,
                "bad component audit")
        require(sha256_file(args.component_results) == audit["results_sha256"],
                "component results hash mismatch")
        arithmetic = json.loads(args.arithmetic_inventory.read_text())
        components = json.loads(args.component_inventory.read_text())
        require(arithmetic.get("schema") == "NONAFFINE_DIAGONAL_1E18_ARITHMETIC_INVENTORY_V1" and
                len(arithmetic.get("cd_shapes", [])) == 39, "bad arithmetic inventory")
        require(components.get("schema") == "CD18_COMPONENT_QUOTIENT_INVENTORY_V2" and
                components.get("class_count") == 381 and
                components.get("arithmetic_inventory_sha256") == sha256_file(args.arithmetic_inventory),
                "bad component inventory")
        rows = [json.loads(line) for line in args.component_results.read_text().splitlines() if line]
        require(len(rows) == 381 and len({row["task_id"] for row in rows}) == 381,
                "component result census mismatch")
        by_pair: dict[tuple[int, int], list[dict[str, Any]]] = {}
        for row in rows:
            by_pair.setdefault((int(row["simple_id"]), int(row["k"])), []).append(row)

        groups = top_groups()
        top_records = []
        top_data: dict[tuple[int, str], tuple[dict[int, int], int]] = {}
        for degree, entries in groups.items():
            for name, group, order in entries:
                coefficients = distinguishing_coefficients(group)
                require(coefficients, f"no distinguishing partition: {name}")
                distinguishing_number = min(coefficients)
                top_data[(degree, name)] = (coefficients, distinguishing_number)
                top_records.append({"m": degree, "name": name, "order": order,
                                    "distinguishing_number": distinguishing_number,
                                    "partition_coefficients": {str(key): value
                                                               for key, value in sorted(coefficients.items())}})

        results = []
        for shape in arithmetic["cd_shapes"]:
            sid, k, m = (int(shape[key]) for key in ("simple_id", "k", "m"))
            degree = int(shape["component_degree"])
            compound_degree = int(shape["degree"])
            require(compound_degree == degree**m and (sid, k) in by_pair,
                    "shape/component mismatch")
            for component in by_pair[(sid, k)]:
                require(int(component["component_degree"]) == degree, "component degree drift")
                bound = Fraction(int(component["bound_num"]), int(component["bound_den"]))
                h_order = int(component["H_order"])
                regular_lower = (0 if bound >= 1 else
                                 ceil_fraction(Fraction(degree, h_order) * (1 - bound)))
                for top_name, _, top_order in groups[m]:
                    coefficients, distinguishing_number = top_data[(m, top_name)]
                    colourings = distinguishing_colourings(coefficients, regular_lower)
                    density = Fraction(colourings * h_order**m, compound_degree)
                    common_regular_colours = max(
                        0,
                        ceil_fraction(Fraction(2 * regular_lower * h_order - degree,
                                               h_order)),
                    )
                    closed_half = density > Fraction(1, 2)
                    closed_distinct = common_regular_colours >= m
                    closed = closed_half or closed_distinct
                    results.append({
                        "overgroup_id": f"CD18_sid{sid}_k{k}_m{m}_case{component['case']}_{top_name}",
                        "simple_id": sid,
                        "simple_name": component["simple_name"],
                        "simple_order": int(component["simple_order"]),
                        "k": k,
                        "m": m,
                        "component_degree": degree,
                        "compound_degree": compound_degree,
                        "component_case": int(component["case"]),
                        "component_B_order": int(component["B_order"]),
                        "component_H_order": h_order,
                        "component_qhat_num": bound.numerator,
                        "component_qhat_den": bound.denominator,
                        "regular_orbits_lower": regular_lower,
                        "top_name": top_name,
                        "top_order": top_order,
                        "distinguishing_number": distinguishing_number,
                        "distinguishing_colourings_lower": colourings,
                        "density_num": density.numerator,
                        "density_den": density.denominator,
                        "common_regular_colours_lower": common_regular_colours,
                        "closed_by_half_density": closed_half,
                        "closed_by_distinct_regular_colours": closed_distinct,
                        "closed": closed,
                    })

        expected = sum(len(by_pair[(int(shape["simple_id"]), int(shape["k"]))]) *
                       len(groups[int(shape["m"])]) for shape in arithmetic["cd_shapes"])
        require(len(results) == expected and len({row["overgroup_id"] for row in results}) == expected,
                "wreath overgroup census mismatch")
        residual = sorted(row["overgroup_id"] for row in results
                          if not row["closed"])
        output = fresh(args.output_dir)
        result_path = output / "CD18_WREATH_RESULTS.jsonl"
        create(result_path, b"".join(canonical(row) for row in sorted(
            results, key=lambda row: row["overgroup_id"])))
        create(output / "CD18_TOP_GROUPS.json", canonical({
            "schema": "CD18_TRANSITIVE_TOP_GROUPS_V1", "complete": True,
            "class_count": 13, "records": top_records,
        }))
        value = {
            "schema": "CD18_WREATH_ANALYSIS_V2",
            "complete": len(residual) == 0,
            "overgroup_count": expected,
            "closed_count": expected - len(residual),
            "residual": residual,
            "component_audit_sha256": sha256_file(args.component_audit),
            "component_results_sha256": sha256_file(args.component_results),
            "arithmetic_inventory_sha256": sha256_file(args.arithmetic_inventory),
            "component_inventory_sha256": sha256_file(args.component_inventory),
            "results_sha256": sha256_file(result_path),
            "proof_rules": [
                "exact_distinguishing_colouring_count_and_strict_half_neighbour_density",
                "regular_colour_intersection_bound_plus_distinct_representatives",
            ],
        }
        create(output / "CD18_WREATH_ANALYSIS.json", canonical(value))
        create(output / "CD18_WREATH_ANALYSIS.sha256",
               f"{sha256_file(output / 'CD18_WREATH_ANALYSIS.json')}  CD18_WREATH_ANALYSIS.json\n".encode())
        print(f"CD18_WREATH_ANALYSIS_OK closed={expected-len(residual)} residual={len(residual)}")
        return 0
    except (OSError, ValueError, json.JSONDecodeError, AnalysisError) as error:
        print(f"CD18_WREATH_ANALYSIS_ERROR {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
