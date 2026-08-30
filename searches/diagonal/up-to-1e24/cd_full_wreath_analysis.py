#!/usr/bin/env python3
"""Close CD24 actions through the containing full symmetric wreath product."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any


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


def create(path: Path, content: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o400)
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


def ceil_fraction(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def falling(value: int, length: int) -> int:
    return 0 if value < length else math.prod(range(value - length + 1, value + 1))


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
        require(audit.get("schema") == "CD24_COMPONENT_PAIR_FPR_AUDIT_V1" and
                audit.get("complete") is True and audit.get("task_count") == 836,
                "bad component audit")
        require(sha256_file(args.component_results) == audit["results_sha256"],
                "component results hash mismatch")
        arithmetic = json.loads(args.arithmetic_inventory.read_text())
        components = json.loads(args.component_inventory.read_text())
        require(arithmetic.get("schema") ==
                "NONAFFINE_DIAGONAL_1E24_ARITHMETIC_INVENTORY_V1" and
                len(arithmetic.get("cd_shapes", [])) == 69, "bad arithmetic inventory")
        require(components.get("schema") == "CD24_COMPONENT_QUOTIENT_INVENTORY_V1" and
                components.get("class_count") == 836 and
                components.get("arithmetic_inventory_sha256") ==
                sha256_file(args.arithmetic_inventory), "bad component inventory")
        rows = [json.loads(line) for line in args.component_results.read_text().splitlines()
                if line]
        require(len(rows) == 836 and len({row["task_id"] for row in rows}) == 836,
                "component result census mismatch")
        by_pair: dict[tuple[int, int], list[dict[str, Any]]] = {}
        for row in rows:
            by_pair.setdefault((int(row["simple_id"]), int(row["k"])), []).append(row)
        expected_pairs = {(int(row["simple_id"]), int(row["k"]))
                          for row in components["component_pairs"]}
        require(set(by_pair) == expected_pairs and len(by_pair) == 69,
                "component pair coverage mismatch")

        results = []
        for shape in arithmetic["cd_shapes"]:
            sid, k, m = (int(shape[key]) for key in ("simple_id", "k", "m"))
            degree = int(shape["component_degree"])
            compound_degree = int(shape["degree"])
            require(compound_degree == degree**m and (sid, k) in by_pair,
                    "shape/component mismatch")
            for component in by_pair[(sid, k)]:
                require(int(component["component_degree"]) == degree,
                        "component degree drift")
                qhat = Fraction(int(component["bound_num"]), int(component["bound_den"]))
                h_order = int(component["H_order"])
                regular_lower = (0 if qhat >= 1 else
                                 ceil_fraction(Fraction(degree, h_order) * (1 - qhat)))
                # A colouring is distinguishing for S_m iff all m colours are distinct.
                colourings = falling(regular_lower, m)
                density = Fraction(colourings * h_order**m, compound_degree)
                common_regular_colours = max(
                    0, ceil_fraction(Fraction(2 * regular_lower * h_order - degree,
                                              h_order)))
                half = density > Fraction(1, 2)
                distinct = common_regular_colours >= m
                results.append({
                    "overgroup_id": f"CD24_sid{sid}_k{k}_m{m}_case{component['case']}_S{m}",
                    "simple_id": sid, "simple_name": component["simple_name"],
                    "simple_order": int(component["simple_order"]), "k": k, "m": m,
                    "component_degree": degree, "compound_degree": compound_degree,
                    "component_case": int(component["case"]),
                    "component_B_order": int(component["B_order"]),
                    "component_H_order": h_order,
                    "component_qhat_num": qhat.numerator,
                    "component_qhat_den": qhat.denominator,
                    "regular_orbits_lower": regular_lower,
                    "top_name": f"S{m}", "top_order": math.factorial(m),
                    "distinguishing_colourings_lower": colourings,
                    "density_num": density.numerator, "density_den": density.denominator,
                    "common_regular_colours_lower": common_regular_colours,
                    "closed_by_full_wreath_half_density": half,
                    "closed_by_distinct_regular_colours": distinct,
                    "closed": half or distinct,
                })
        expected = sum(len(by_pair[(int(shape["simple_id"]), int(shape["k"]))])
                       for shape in arithmetic["cd_shapes"])
        require(expected == 836 and len(results) == expected and
                len({row["overgroup_id"] for row in results}) == expected,
                "full-wreath overgroup census mismatch")
        residual = sorted(row["overgroup_id"] for row in results if not row["closed"])
        output = fresh(args.output_dir)
        result_path = output / "CD24_FULL_WREATH_RESULTS.jsonl"
        create(result_path, b"".join(canonical(row) for row in sorted(
            results, key=lambda row: row["overgroup_id"])))
        value = {
            "schema": "CD24_FULL_WREATH_ANALYSIS_V1",
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
                "full_symmetric_top_exact_distinguishing_colourings_and_strict_half_density",
                "regular_colour_intersection_bound_plus_distinct_representatives",
                "closure_in_full_wreath_implies_closure_for_every_contained_cd_group",
            ],
        }
        audit_path = output / "CD24_FULL_WREATH_ANALYSIS.json"
        create(audit_path, canonical(value))
        create(output / "CD24_FULL_WREATH_ANALYSIS.sha256",
               f"{sha256_file(audit_path)}  CD24_FULL_WREATH_ANALYSIS.json\n".encode())
        print(f"CD24_FULL_WREATH_ANALYSIS_OK closed={expected-len(residual)} "
              f"residual={len(residual)}")
        return 0
    except (OSError, ValueError, json.JSONDecodeError, AnalysisError) as error:
        print(f"CD24_FULL_WREATH_ANALYSIS_ERROR {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
