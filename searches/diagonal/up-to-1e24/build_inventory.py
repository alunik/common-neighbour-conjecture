#!/usr/bin/env python3
"""Build and audit the canonical 10^18 < n <= 10^24 SD shape inventory."""

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
SIMPLE_LIMIT = 10**12


def require(value: bool, message: str) -> None:
    if not value:
        raise RuntimeError(message)


def parse(line: str) -> tuple[str, dict[str, str]]:
    parts = line.rstrip("\n").split("|")
    fields: dict[str, str] = {}
    for item in parts[1:]:
        key, value = item.split("=", 1)
        fields[key] = value
    return parts[0], fields


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


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--magma", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output_dir.exists() and not args.output_dir.is_symlink(),
            "output must be fresh")
    completed = subprocess.run([str(args.magma), str(args.source)], check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True)
    require(completed.stderr == "", f"unexpected Magma stderr: {completed.stderr}")
    simples: list[dict[str, Any]] = []
    sd: list[dict[str, Any]] = []
    summary: dict[str, int] | None = None
    cutoff: dict[str, int] | None = None
    for line in completed.stdout.splitlines():
        if "|" not in line:
            continue
        kind, row = parse(line)
        if kind == "SIMPLE":
            simples.append({"simple_id": int(row["sid"]), "name": row["name"],
                            "order": int(row["order"])})
        elif kind == "SD_SHAPE":
            sd.append({"simple_id": int(row["sid"]), "name": row["name"],
                       "order": int(row["order"]), "k": int(row["k"]),
                       "degree": int(row["degree"])})
        elif kind == "CATALOGUE_CUTOFF":
            cutoff = {key: int(value) for key, value in row.items()}
        elif kind == "SHAPE_SUMMARY":
            summary = {key: int(value) for key, value in row.items()}
    require(summary is not None and cutoff is not None, "missing inventory terminal records")
    require(len(simples) == summary["simple_groups"] and
            [row["simple_id"] for row in simples] == list(range(1, len(simples) + 1)) and
            simples[-1]["order"] <= SIMPLE_LIMIT < cutoff["next_order"] and
            cutoff["next_sid"] == len(simples) + 1 and cutoff["limit"] == SIMPLE_LIMIT,
            "simple catalogue cutoff mismatch")
    require(len(sd) == summary["sd_shapes"] and
            all(LOWER < int(row["degree"]) <= UPPER for row in sd),
            "shape census/window mismatch")
    require(summary == {"simple_groups": 1650, "sd_shapes": 1556, "max_k": 14},
            "unexpected SD census")
    simple_by_id = {int(row["simple_id"]): row for row in simples}
    expected_sd = set()
    for sid, simple in simple_by_id.items():
        order = int(simple["order"])
        component = order**2
        k = 3
        while component <= UPPER:
            if component > LOWER:
                expected_sd.add((sid, k, component))
            component *= order
            k += 1
    observed_sd = {(int(row["simple_id"]), int(row["k"]), int(row["degree"]))
                   for row in sd}
    require(observed_sd == expected_sd, "independent arithmetic replay mismatch")
    require(summary["max_k"] == max(int(row["k"]) for row in sd),
            "maximum parameter mismatch")
    value = {
        "schema": "DIAGONAL_SD_1E24_ARITHMETIC_INVENTORY_V1",
        "degree_window": {"minimum_exclusive": LOWER, "maximum_inclusive": UPPER},
        "simple_order_limit": SIMPLE_LIMIT,
        "catalogue_cutoff": cutoff,
        "simple_groups": simples,
        "sd_shapes": sorted(sd, key=lambda row: (row["degree"], row["simple_id"], row["k"])),
        "counts": summary,
    }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    transcript = args.output_dir / "MAGMA_TRANSCRIPT.txt"
    inventory = args.output_dir / "ARITHMETIC_INVENTORY.json"
    create(transcript, completed.stdout.encode())
    create(inventory, canonical(value))
    create(args.output_dir / "ARITHMETIC_INVENTORY.sha256",
           f"{sha256_file(inventory)}  ARITHMETIC_INVENTORY.json\n".encode())
    print("DIAGONAL_SD_1E24_INVENTORY_OK "
          f"simple={len(simples)} sd={len(sd)} "
          f"sha256={sha256_file(inventory)}")


if __name__ == "__main__":
    main()
