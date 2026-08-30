#!/usr/bin/env python3
"""Build and cross-audit the CD component quotient inventory."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Any


def parse(line: str) -> tuple[str, dict[str, str]]:
    fields = line.split("|")
    values: dict[str, str] = {}
    for field in fields[1:]:
        key, sep, value = field.partition("=")
        if sep != "=" or not key or key in values:
            raise RuntimeError("bad Magma transcript field")
        values[key] = value
    return fields[0], values


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def create(path: Path, value: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o400)
    try:
        os.write(descriptor, value)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--arithmetic-inventory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    completed = subprocess.run(["magma", "-b", str(args.source.resolve(strict=True))],
                               check=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    if completed.stderr:
        raise RuntimeError(f"unexpected Magma stderr: {completed.stderr}")
    classes: list[dict[str, Any]] = []
    pairs: list[dict[str, Any]] = []
    summary = None
    for line in completed.stdout.splitlines():
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
    if summary != {"pairs": 31, "classes": 381} or len(pairs) != 31 or len(classes) != 381:
        raise RuntimeError("CD component census mismatch")
    arithmetic = json.loads(args.arithmetic_inventory.read_text())
    shape_pairs: dict[tuple[int, int], set[int]] = {}
    for shape in arithmetic["cd_shapes"]:
        shape_pairs.setdefault((int(shape["simple_id"]), int(shape["k"])), set()).add(int(shape["m"]))
    if set(shape_pairs) != {(row["simple_id"], row["k"]) for row in pairs}:
        raise RuntimeError("CD pair coverage mismatch")
    for row in pairs:
        row["compound_m_values"] = sorted(shape_pairs[(row["simple_id"], row["k"])])
        selected = [item for item in classes
                    if (item["simple_id"], item["k"]) == (row["simple_id"], row["k"])]
        if len(selected) != row["class_count"] or [item["case"] for item in selected] != list(
                range(1, row["class_count"] + 1)):
            raise RuntimeError("per-pair class census mismatch")
    result = {"schema": "CD18_COMPONENT_QUOTIENT_INVENTORY_V2", "complete": True,
              "degree_window": {"minimum_exclusive": 100_000_000,
                                "maximum_inclusive": 10**18},
              "pair_count": 31, "class_count": 381,
              "component_pairs": pairs, "component_classes": classes,
              "arithmetic_inventory_sha256": hashlib.sha256(
                  args.arithmetic_inventory.read_bytes()).hexdigest()}
    payload = canonical(result)
    output = args.output.resolve()
    if output.exists() or output.is_symlink():
        raise RuntimeError("output exists")
    output.mkdir(mode=0o700)
    create(output / "MAGMA_TRANSCRIPT.txt", completed.stdout.encode())
    create(output / "CD_COMPONENT_INVENTORY.json", payload)
    create(output / "CD_COMPONENT_INVENTORY.sha256",
           f"{sha256_bytes(payload)}  CD_COMPONENT_INVENTORY.json\n".encode())
    print(f"CD18_COMPONENT_INVENTORY_OK pairs=31 classes=381 sha256={sha256_bytes(payload)}")


if __name__ == "__main__":
    main()
