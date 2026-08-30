#!/usr/bin/env python3
"""Combine completed catalogue-frontier shards into the final census."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, action="append", required=True)
    parser.add_argument("--root-counts", type=Path, required=True)
    parser.add_argument("--dimension-one", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    count_lines = [line.strip() for line in args.root_counts.read_text().splitlines() if line.strip()]
    count_header = count_lines[0].split("\t")
    root_counts = {
        int(row["root_index"]): row
        for row in (
            dict(zip(count_header, line.split("\t"), strict=True))
            for line in count_lines[1:]
        )
    }

    covered_roots: set[int] = set()
    result_rows: list[dict[str, object]] = []
    run_summaries: list[dict[str, object]] = []
    for run_dir_arg in args.run_dir:
        run_dir = run_dir_arg.resolve()
        state_path = run_dir / "state.json"
        nodes_path = run_dir / "nodes.jsonl"
        state = json.loads(state_path.read_text())
        if state.get("status") != "complete" or state.get("queue"):
            raise SystemExit(f"run is not complete: {run_dir}")
        first = int(state["root_first"])
        last = int(state["root_last"])
        roots = set(range(first, last + 1))
        overlap = covered_roots & roots
        if overlap:
            raise SystemExit(f"overlapping root ranges at {min(overlap)}")
        covered_roots |= roots
        rows = [json.loads(line) for line in nodes_path.read_text().splitlines() if line.strip()]
        if len(rows) != int(state["nodes"]):
            raise SystemExit(f"node count mismatch in {run_dir}")
        if int(state["guard_failures"]) != 0:
            raise SystemExit(f"resource guard failures in {run_dir}")
        result_rows.extend(rows)
        run_summaries.append(
            {
                "path": str(run_dir),
                "root_first": first,
                "root_last": last,
                "tested_nodes": len(rows),
                "safe_nodes": int(state["safe_nodes"]),
                "explicit_no_base_two": int(state["no_base_two_nodes"]),
                "diameter_at_least_three": int(state["candidate_nodes"]),
                "state_sha256": sha256(state_path),
                "nodes_sha256": sha256(nodes_path),
                "engine_binary_sha256": sha256(run_dir / "bin" / "affine_saxl_engine"),
                "engine_source_sha256": sha256(run_dir / "source" / "affine_saxl_engine.cpp"),
                "exporter_source_sha256": sha256(
                    run_dir / "source" / "export_irredsol_catalogue_round.g"
                ),
                "controller_source_sha256": sha256(
                    run_dir / "source" / "run_irredsol_catalogue_rounds.py"
                ),
            }
        )

    expected_roots = set(root_counts)
    if covered_roots != expected_roots:
        missing = sorted(expected_roots - covered_roots)
        extra = sorted(covered_roots - expected_roots)
        raise SystemExit(f"root coverage mismatch: missing={missing[:10]} extra={extra[:10]}")

    labels = [str(row["label"]) for row in result_rows]
    if len(labels) != len(set(labels)):
        raise SystemExit("duplicate tested catalogue labels across runs")
    no_base = [row for row in result_rows if row.get("status") == "no_regular_vector"]
    candidates = [row for row in result_rows if row.get("status") == "diameter_at_least_3"]
    unexpected = [
        row for row in result_rows
        if row.get("status") not in {
            "no_regular_vector", "diameter_at_least_3", "density_certificate",
            "fixed_space_union_bound", "two_sum_certificate", "exact",
        }
    ]
    if unexpected:
        raise SystemExit(f"unexpected result statuses: {unexpected[:3]}")

    n_ge_2_actions = sum(int(row["actions"]) for row in root_counts.values())
    n_ge_2_eligible = sum(int(row["order_eligible"]) for row in root_counts.values())
    n_ge_2_order_obstructed = sum(int(row["order_obstructed"]) for row in root_counts.values())
    n_ge_2_base_two = n_ge_2_eligible - len(no_base)
    n_ge_2_diameter_at_most_two = n_ge_2_base_two - len(candidates)

    dimension_one = json.loads(args.dimension_one.read_text())
    all_actions = int(dimension_one["actions"]) + n_ge_2_actions
    all_base_one = int(dimension_one["base_size_1"])
    all_base_two = int(dimension_one["base_size_2"]) + n_ge_2_base_two
    all_base_gt_two = n_ge_2_order_obstructed + len(no_base)
    if all_base_one + all_base_two + all_base_gt_two != all_actions:
        raise SystemExit("base-size totals do not partition all actions")

    summary = {
        "schema": "IRREDSOL_AFFINE_CENSUS_SUMMARY_V1",
        "degree_lower": 8192,
        "degree_upper": 2**24 - 1,
        "runs": run_summaries,
        "root_count": len(root_counts),
        "tested_frontier_nodes": len(result_rows),
        "dimension_at_least_two": {
            "actions": n_ge_2_actions,
            "order_eligible": n_ge_2_eligible,
            "order_obstructed_base_gt_two": n_ge_2_order_obstructed,
            "exact_no_regular_vector": len(no_base),
            "base_size_two": n_ge_2_base_two,
            "diameter_at_most_two": n_ge_2_diameter_at_most_two,
            "diameter_at_least_three": len(candidates),
        },
        "dimension_one": dimension_one,
        "all_dimensions": {
            "actions": all_actions,
            "base_size_one": all_base_one,
            "base_size_two": all_base_two,
            "base_size_greater_than_two": all_base_gt_two,
            "diameter_at_most_two_among_base_two": (
                int(dimension_one["diameter_at_most_2"]) + n_ge_2_diameter_at_most_two
            ),
            "diameter_at_least_three_among_base_two": len(candidates),
        },
        "diameter_at_least_three_actions": candidates,
        "root_counts_sha256": sha256(args.root_counts),
        "dimension_one_sha256": sha256(args.dimension_one),
    }
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["all_dimensions"], sort_keys=True))


if __name__ == "__main__":
    main()
