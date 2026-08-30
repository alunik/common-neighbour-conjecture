#!/usr/bin/env python3
"""Extract immutable per-case Magma descriptors from an exact frontier transcript."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any


class DataError(RuntimeError):
    pass


def canonical(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def create(path: Path, content: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o400)
    try:
        os.write(descriptor, content)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fields(line: str) -> dict[str, str]:
    answer = {}
    for part in line.rstrip("\n").split("|")[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            answer[key] = value
    return answer


def numbers(text: str) -> list[int]:
    return [int(value) for value in re.findall(r"[0-9]+", text)]


def magma_sequence(values: list[Any], depth: int = 0) -> str:
    if not values:
        return "[]"
    if isinstance(values[0], list):
        indent = "    " * depth
        child = ",\n".join("    " * (depth + 1) + magma_sequence(value, depth + 1)
                            for value in values)
        return "[\n" + child + "\n" + indent + "]"
    return "[" + ",".join(str(value) for value in values) + "]"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frontier", type=Path, required=True)
    parser.add_argument("--mode", choices=("full", "s3"), required=True)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--case", type=int, action="append")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    frontier = args.frontier.resolve(strict=True)
    cases: dict[int, dict[str, Any]] = {}
    terminal = None
    for line in frontier.read_text().splitlines():
        if line.startswith("FRONTIER_BASE2|"):
            row = fields(line)
            if row.get("gt_half") == "false":
                number = int(row["node"])
                cases[number] = {"case": number, "m": 6, "order": int(row["order"]),
                                 "regular": int(row["regular"]), "tops": [],
                                 "outers": [], "perms": []}
        elif line.startswith("FRONTIER_GENERATOR|"):
            row = fields(line)
            number = int(row["node"])
            if number not in cases:
                continue
            top_start = line.index("|top=") + len("|top=")
            components_start = line.index("|components=")
            top = numbers(line[top_start:components_start])
            component_text = line[components_start + len("|components=") + 1:-1]
            components = component_text.split(";")
            outer_bits, permutations = [], []
            for component in components:
                bit_text, permutation_text = component.split(":", 1)
                outer_bits.append(int(bit_text))
                permutations.append(numbers(permutation_text))
            if len(top) != 6 or len(outer_bits) != 6 or any(len(p) != 3 for p in permutations):
                raise DataError(f"bad generator for case {number}")
            cases[number]["tops"].append(top)
            cases[number]["outers"].append(outer_bits)
            cases[number]["perms"].append(permutations)
        elif line.startswith("FRONTIER_COMPLETE|"):
            terminal = fields(line)
    if not args.allow_incomplete:
        if terminal is None or terminal.get("component") != args.mode or terminal.get("m") != "6":
            raise DataError("frontier is not terminal for the requested mode")
        if int(terminal["base2_frontier"]) < len(cases):
            raise DataError("frontier residual census mismatch")
    selected = sorted(set(args.case or cases))
    if not selected or any(number not in cases for number in selected):
        raise DataError("requested case is absent")
    output = args.output_dir.resolve()
    if output.exists() or output.is_symlink():
        raise DataError("output exists")
    output.mkdir(parents=True, mode=0o700)
    records = []
    for number in selected:
        item = cases[number]
        count = len(item["tops"])
        if count == 0 or len(item["outers"]) != count or len(item["perms"]) != count:
            raise DataError(f"incomplete generators for case {number}")
        text = (
            f'CaseMode := "{args.mode}";\n'
            f'CaseNumber := {number};\n'
            f'CaseExpectedOrder := {item["order"]};\n'
            f'CaseExpectedRegular := {item["regular"]};\n'
            f'CaseTops := {magma_sequence(item["tops"])};\n'
            f'CaseOuters := {magma_sequence(item["outers"])};\n'
            f'CasePerms := {magma_sequence(item["perms"])};\n')
        path = output / f"case_{number:04d}.m"
        create(path, text.encode())
        records.append({"case": number, "order": item["order"],
                        "regular": item["regular"], "generators": count,
                        "file": path.name, "sha256": sha256_file(path)})
    manifest = {"schema": "DIAGONAL24_CD_FRONTIER_CASE_DATA_V1", "complete": True,
                "mode": args.mode, "m": 6, "frontier": str(frontier),
                "frontier_sha256": sha256_file(frontier), "case_count": len(records),
                "frontier_terminal": terminal, "cases": records}
    manifest_path = output / "CASE_DATA_MANIFEST.json"
    create(manifest_path, canonical(manifest))
    create(output / "CASE_DATA_MANIFEST.sha256",
           f"{sha256_file(manifest_path)}  CASE_DATA_MANIFEST.json\n".encode())
    print(f"DIAGONAL24_CD_CASE_DATA_OK mode={args.mode} cases={len(records)} "
          f"sha256={sha256_file(manifest_path)}")


if __name__ == "__main__":
    main()
