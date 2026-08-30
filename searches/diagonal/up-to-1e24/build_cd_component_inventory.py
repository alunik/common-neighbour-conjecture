#!/usr/bin/env python3
"""Build and independently cross-audit the CD component quotient inventory."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Any


LOWER = 10**18
UPPER = 10**24


def require(value: bool, message: str) -> None:
    if not value:
        raise RuntimeError(message)


def parse(line: str) -> tuple[str, dict[str, str]]:
    fields = line.split("|")
    values: dict[str, str] = {}
    for field in fields[1:]:
        key, sep, value = field.partition("=")
        require(sep == "=" and bool(key) and key not in values,
                "bad Magma transcript field")
        values[key] = value
    return fields[0], values


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def create(path: Path, value: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o400)
    try:
        os.write(descriptor, value)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--magma", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--arithmetic-inventory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(),
            "output must be fresh")
    completed = subprocess.run([str(args.magma), "-b", str(args.source.resolve(strict=True))],
                               check=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    require(completed.stderr == "", f"unexpected Magma stderr: {completed.stderr}")
    classes: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    summary: dict[str, int] | None = None
    for line in completed.stdout.splitlines():
        if "|" not in line:
            continue
        kind, row = parse(line)
        if kind == "CD_COMPONENT_CLASS":
            classes.append({"simple_id": int(row["sid"]), "name": row["name"],
                            "order": int(row["order"]), "k": int(row["k"]),
                            "component_degree": int(row["component_degree"]),
                            "case": int(row["case"]), "B_order": int(row["B"]),
                            "top_order": int(row["top"]),
                            "outer_image_order": int(row["outer_image"]),
                            "top_kernel_order": int(row["top_kernel"]),
                            "outer_order": int(row["Out"])})
        elif kind == "CD_COMPONENT_PAIR":
            pairs.append({"simple_id": int(row["sid"]), "name": row["name"],
                          "order": int(row["order"]), "k": int(row["k"]),
                          "component_degree": int(row["component_degree"]),
                          "class_count": int(row["classes"]),
                          "outer_order": int(row["Out"])})
        elif kind == "CD_COMPONENT_SUMMARY":
            summary = {key: int(value) for key, value in row.items()}
    require(summary is not None and summary["pairs"] == len(pairs) and
            summary["classes"] == len(classes), "CD component census mismatch")
    arithmetic = json.loads(args.arithmetic_inventory.read_text())
    require(arithmetic["schema"] == "NONAFFINE_DIAGONAL_1E24_ARITHMETIC_INVENTORY_V1" and
            arithmetic["degree_window"] == {
                "minimum_exclusive": LOWER, "maximum_inclusive": UPPER},
            "arithmetic inventory contract mismatch")
    shape_pairs: dict[tuple[int, int], set[int]] = {}
    simple_orders = {int(row["simple_id"]): int(row["order"])
                     for row in arithmetic["simple_groups"]}
    for shape in arithmetic["cd_shapes"]:
        key = (int(shape["simple_id"]), int(shape["k"]))
        shape_pairs.setdefault(key, set()).add(int(shape["m"]))
        require(int(shape["component_degree"]) == simple_orders[key[0]] ** (key[1] - 1),
                "component degree mismatch")
    observed_pairs = {(row["simple_id"], row["k"]) for row in pairs}
    require(len(pairs) == len(observed_pairs) == len(shape_pairs) and
            observed_pairs == set(shape_pairs), "CD pair coverage mismatch")
    for row in pairs:
        key = (row["simple_id"], row["k"])
        require(row["order"] == simple_orders[key[0]] and
                row["component_degree"] == row["order"] ** (row["k"] - 1),
                "pair arithmetic mismatch")
        row["compound_m_values"] = sorted(shape_pairs[key])
        selected = [item for item in classes
                    if (item["simple_id"], item["k"]) == key]
        require(len(selected) == row["class_count"] and
                [item["case"] for item in selected] == list(
                    range(1, row["class_count"] + 1)),
                "per-pair class census mismatch")
        require(all(item["order"] == row["order"] and
                    item["component_degree"] == row["component_degree"] and
                    item["outer_order"] == row["outer_order"] for item in selected),
                "per-class pair binding mismatch")
    result = {
        "schema": "CD24_COMPONENT_QUOTIENT_INVENTORY_V1",
        "complete": True,
        "degree_window": {"minimum_exclusive": LOWER,
                          "maximum_inclusive": UPPER},
        "pair_count": len(pairs),
        "class_count": len(classes),
        "component_pairs": pairs,
        "component_classes": classes,
        "arithmetic_inventory_sha256": hashlib.sha256(
            args.arithmetic_inventory.read_bytes()).hexdigest(),
        "magma_source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
    }
    payload = canonical(result)
    args.output.mkdir(mode=0o700)
    create(args.output / "MAGMA_TRANSCRIPT.txt", completed.stdout.encode())
    create(args.output / "CD_COMPONENT_INVENTORY.json", payload)
    create(args.output / "CD_COMPONENT_INVENTORY.sha256",
           f"{sha256_bytes(payload)}  CD_COMPONENT_INVENTORY.json\n".encode())
    print("CD24_COMPONENT_INVENTORY_OK "
          f"pairs={len(pairs)} classes={len(classes)} sha256={sha256_bytes(payload)}")


if __name__ == "__main__":
    main()
