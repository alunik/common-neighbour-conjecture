#!/usr/bin/env python3
"""Build a canonical JSON inventory from enumerate_shapes.m output."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def parse_fields(line: str) -> tuple[str, dict[str, str]]:
    fields = line.rstrip("\n").split("|")
    kind = fields[0]
    values: dict[str, str] = {}
    for field in fields[1:]:
        key, value = field.split("=", 1)
        values[key] = value
    return kind, values


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("transcript", type=Path, nargs="?")
    parser.add_argument("--magma-source", type=Path)
    args = parser.parse_args()

    if (args.transcript is None) == (args.magma_source is None):
        parser.error("provide exactly one transcript or --magma-source")
    if args.magma_source is not None:
        completed = subprocess.run(
            ["magma", str(args.magma_source)], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if completed.stderr:
            raise RuntimeError(f"unexpected Magma stderr: {completed.stderr}")
        transcript_text = completed.stdout
    else:
        transcript_text = args.transcript.read_text()

    simples: list[dict[str, object]] = []
    sd: list[dict[str, object]] = []
    cd: list[dict[str, object]] = []
    summary: dict[str, int] | None = None

    for line in transcript_text.splitlines():
        kind, row = parse_fields(line)
        if kind == "SIMPLE":
            simples.append({"simple_id": int(row["sid"]), "name": row["name"],
                            "order": int(row["order"])})
        elif kind == "SD_SHAPE":
            sd.append({"simple_id": int(row["sid"]), "name": row["name"],
                       "order": int(row["order"]), "k": int(row["k"]),
                       "degree": int(row["degree"])})
        elif kind == "CD_SHAPE":
            cd.append({"simple_id": int(row["sid"]), "name": row["name"],
                       "order": int(row["order"]), "k": int(row["k"]),
                       "component_degree": int(row["component_degree"]),
                       "m": int(row["m"]), "degree": int(row["degree"])})
        elif kind == "SHAPE_SUMMARY":
            summary = {key: int(value) for key, value in row.items()}

    assert summary == {"simple_groups": 277, "sd_shapes": 357,
                       "cd_shapes": 39, "max_k": 11, "max_m": 5}
    assert len(simples) == 277 and [r["simple_id"] for r in simples] == list(range(1, 278))
    assert len(sd) == 357 and len(cd) == 39
    assert all(100_000_000 < int(r["degree"]) <= 10**18 for r in sd + cd)
    assert max(int(r["k"]) for r in sd) == 11
    assert max(int(r["m"]) for r in cd) == 5

    inventory = {
        "schema": "NONAFFINE_DIAGONAL_1E18_ARITHMETIC_INVENTORY_V1",
        "degree_window": {"minimum_exclusive": 100_000_000,
                          "maximum_inclusive": 10**18},
        "simple_groups": simples,
        "sd_shapes": sorted(sd, key=lambda r: (int(r["degree"]), int(r["simple_id"]), int(r["k"]))),
        "cd_shapes": sorted(cd, key=lambda r: (int(r["degree"]), int(r["simple_id"]),
                                                int(r["k"]), int(r["m"]))),
        "counts": summary,
        "routing": {
            "sd_non_ak_sk": "huang_thesis_theorem_5_6",
            "sd_ak_sk": "full_normalizer_fpr_then_residual",
            "cd": "component_full_wreath_density_then_exact_rainbow_residual",
        },
    }
    payload = canonical(inventory)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "MAGMA_TRANSCRIPT.txt").write_text(transcript_text)
    inventory_path = args.output / "ARITHMETIC_INVENTORY.json"
    inventory_path.write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    (args.output / "ARITHMETIC_INVENTORY.sha256").write_text(
        f"{digest}  ARITHMETIC_INVENTORY.json\n")
    print(f"DIAGONAL_1E18_INVENTORY_OK sd=357 cd=39 sha256={digest}")


if __name__ == "__main__":
    main()
