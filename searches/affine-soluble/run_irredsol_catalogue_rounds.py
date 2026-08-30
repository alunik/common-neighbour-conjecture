#!/usr/bin/env python3
"""Layer-wise, resumable IRREDSOL search for Saxl diameter at least three."""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ENGINE_SOURCE = ROOT / "searches" / "engines" / "affine_saxl_engine.cpp"
EXPORTER = HERE / "export_irredsol_catalogue_round.g"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def recover_jsonl(path: Path, committed_lines: int) -> None:
    """Discard output appended by an interrupted, uncommitted batch."""
    lines = path.read_text(encoding="utf-8").splitlines() if path.is_file() else []
    if len(lines) < committed_lines:
        raise RuntimeError(
            f"{path} has {len(lines)} lines, below committed count {committed_lines}"
        )
    if len(lines) > committed_lines:
        path.write_text(
            "\n".join(lines[:committed_lines]) + ("\n" if committed_lines else ""),
            encoding="utf-8",
        )


def family_key(task: dict[str, object]) -> str:
    return f"{task['n']}:{task['p']}:{task['d']}:{task['guardian']}"


def task_key(task: dict[str, object]) -> str:
    return f"{family_key(task)}:{task['round']}"


def task_from_label(label: str) -> tuple[dict[str, object], int]:
    fields = label.split("_")
    if len(fields) != 6 or fields[0] != "ICR":
        raise RuntimeError(f"cannot decode catalogue-round label: {label}")
    task = {
        "n": int(fields[1]), "p": int(fields[2]), "d": int(fields[3]),
        "guardian": int(fields[4]),
    }
    return task, int(fields[5])


def gap_record(task: dict[str, object], tested: list[int], safe: list[int]) -> str:
    tested_text = ",".join(str(k) for k in tested)
    safe_text = ",".join(str(k) for k in safe)
    return (
        f"rec(n:={task['n']},p:={task['p']},d:={task['d']},"
        f"guardian:={task['guardian']},tested:=[{tested_text}],safe:=[{safe_text}])"
    )


def load_roots(path: Path, first: int, last: int | None) -> list[dict[str, object]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    header = lines[0].split("\t")
    if header != ["root_index", "degree", "n", "p", "d", "guardian"]:
        raise RuntimeError(f"unexpected root manifest header: {header}")
    roots: list[dict[str, object]] = []
    for line in lines[1:]:
        values = dict(zip(header, line.split("\t"), strict=True))
        index = int(values["root_index"])
        if index < first or (last is not None and index > last):
            continue
        roots.append({
            "n": int(values["n"]), "p": int(values["p"]),
            "d": int(values["d"]), "guardian": int(values["guardian"]),
            "round": 0,
        })
    return roots


def action_blocks(path: Path) -> list[list[str]]:
    lines = [line.rstrip("\r") for line in path.read_text().splitlines() if line.strip()]
    if not lines or lines[0].strip() != "AFFINE_SAXL_V1":
        raise RuntimeError(f"bad affine input header in {path}")
    blocks: list[list[str]] = []
    current: list[str] | None = None
    for line in lines[1:]:
        if line.strip() == "action":
            if current is not None:
                raise RuntimeError(f"nested action block in {path}")
            current = [line]
        elif current is not None:
            current.append(line)
            if line.strip() == "end":
                blocks.append(current)
                current = None
        else:
            raise RuntimeError(f"content outside action block in {path}: {line!r}")
    if current is not None:
        raise RuntimeError(f"unterminated action block in {path}")
    return blocks


def run_engine_parallel(
    engine: Path,
    input_path: Path,
    result_path: Path,
    engine_log: Path,
    jobs: int,
    pair_budget: int,
) -> int:
    blocks = action_blocks(input_path)
    workers = min(jobs, len(blocks))
    if workers == 0:
        return 0
    # More chunks than workers prevents one unusually expensive contiguous
    # slice from leaving the other workers idle near the end of a batch.
    chunk_count = min(len(blocks), 4 * workers)
    parts: list[tuple[Path, Path, Path]] = []
    for part in range(chunk_count):
        first = len(blocks) * part // chunk_count
        last = len(blocks) * (part + 1) // chunk_count
        part_input = input_path.with_name(f"{input_path.stem}.part_{part + 1:02d}.asx")
        part_result = result_path.with_name(f"{result_path.stem}.part_{part + 1:02d}.jsonl")
        part_log = engine_log.with_name(f"{engine_log.stem}.part_{part + 1:02d}.log")
        text = "AFFINE_SAXL_V1\n" + "\n".join(
            "\n".join(block) for block in blocks[first:last]
        ) + "\n"
        part_input.write_text(text)
        parts.append((part_input, part_result, part_log))

    def run_part(paths: tuple[Path, Path, Path]) -> tuple[Path, int]:
        part_input, part_result, part_log = paths
        with part_result.open("w") as output, part_log.open("w") as diagnostics:
            completed = subprocess.run(
                [str(engine), str(part_input), "--diameter-ge3-only",
                 f"--pair-budget={pair_budget}", "--max-vectors=16777215"],
                cwd=ROOT, stdin=subprocess.DEVNULL, text=True,
                stdout=output, stderr=diagnostics,
            )
        return part_result, completed.returncode

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        completed_parts = list(executor.map(run_part, parts))
    failures = [(path, code) for path, code in completed_parts if code != 0]
    if failures:
        raise RuntimeError(f"parallel engine failures: {failures}")
    with result_path.open("w") as destination:
        for part_result, _ in completed_parts:
            destination.write(part_result.read_text())
    engine_log.write_text(
        f"parallel_engine_workers={workers}\nparallel_engine_parts={chunk_count}\n" +
        "\n".join(str(part_log) for _, _, part_log in parts) + "\n"
    )
    return workers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("roots", type=Path)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--root-first", type=int, default=1)
    parser.add_argument("--root-last", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument("--pair-budget", type=int, default=200_000_000)
    parser.add_argument("--engine-jobs", type=int, default=1)
    parser.add_argument("--max-batches", type=int)
    args = parser.parse_args()
    if (args.root_last < args.root_first or args.batch_size <= 0 or
            args.pair_budget < 0 or args.engine_jobs <= 0):
        raise SystemExit("invalid batch size, pair budget, or engine job count")

    run_dir = args.run_dir.resolve()
    task_dir, input_dir = run_dir / "tasks", run_dir / "inputs"
    result_dir, log_dir = run_dir / "results", run_dir / "logs"
    source_dir, bin_dir = run_dir / "source", run_dir / "bin"
    for directory in (task_dir, input_dir, result_dir, log_dir, source_dir, bin_dir):
        directory.mkdir(parents=True, exist_ok=True)

    state_path = run_dir / "state.json"
    snapshots = [
        (ENGINE_SOURCE, source_dir / ENGINE_SOURCE.name),
        (EXPORTER, source_dir / EXPORTER.name),
        (Path(__file__), source_dir / Path(__file__).name),
        (args.roots, source_dir / args.roots.name),
    ]
    if state_path.is_file():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        for current, snapshot in snapshots:
            if not snapshot.is_file() or sha256(current) != sha256(snapshot):
                raise RuntimeError(
                    f"source changed since this run started: {current}; "
                    "start a new run directory"
                )
        expected_configuration = {
            "root_first": args.root_first, "root_last": args.root_last,
            "batch_size": args.batch_size, "pair_budget": args.pair_budget,
            "engine_jobs": args.engine_jobs,
        }
        for key, value in expected_configuration.items():
            if state.get(key) != value:
                raise RuntimeError(
                    f"run configuration changed for {key}: "
                    f"stored={state.get(key)!r} requested={value!r}"
                )
    else:
        for current, snapshot in snapshots:
            shutil.copy2(current, snapshot)
        roots = load_roots(args.roots, args.root_first, args.root_last)
        state = {
            "schema": "IRREDSOL_CATALOGUE_ROUNDS_STATE_V1",
            "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "root_first": args.root_first, "root_last": args.root_last,
            "root_count": len(roots), "queue": roots,
            "batch_size": args.batch_size, "pair_budget": args.pair_budget,
            "engine_jobs": args.engine_jobs,
            "processed_task_keys": [], "seen_node_labels": [],
            "tested_indices": {}, "safe_indices": {},
            "batch": 0, "nodes": 0, "safe_nodes": 0,
            "no_base_two_nodes": 0, "candidate_nodes": 0,
            "empty_rounds": 0, "guard_failures": 0,
        }
        atomic_json(state_path, state)

    engine = bin_dir / "affine_saxl_engine"
    run_engine_source = source_dir / ENGINE_SOURCE.name
    run_exporter = source_dir / EXPORTER.name
    subprocess.run(
        ["clang++", "-O3", "-march=native", "-DNDEBUG", "-std=c++17",
         str(run_engine_source), "-o", str(engine)], cwd=ROOT, check=True,
    )

    processed = set(state["processed_task_keys"])
    seen_labels = set(state["seen_node_labels"])
    tested_indices: dict[str, list[int]] = state["tested_indices"]
    safe_indices: dict[str, list[int]] = state["safe_indices"]
    master_results = run_dir / "nodes.jsonl"
    candidates_path = run_dir / "candidates.jsonl"
    shard_manifest = run_dir / "shards.jsonl"
    recover_jsonl(master_results, int(state["nodes"]))
    recover_jsonl(candidates_path, int(state["candidate_nodes"]))
    recover_jsonl(shard_manifest, int(state["batch"]))
    batches_this_call = 0

    while state["queue"]:
        if args.max_batches is not None and batches_this_call >= args.max_batches:
            break
        batch_tasks: list[dict[str, object]] = []
        while state["queue"] and len(batch_tasks) < args.batch_size:
            task = state["queue"].pop(0)
            key = task_key(task)
            if key in processed:
                continue
            processed.add(key)
            batch_tasks.append(task)
        if not batch_tasks:
            continue

        state["batch"] += 1
        batches_this_call += 1
        stem = f"batch_{state['batch']:06d}"
        task_path, input_path = task_dir / f"{stem}.g", input_dir / f"{stem}.asx"
        result_path = result_dir / f"{stem}.jsonl"
        gap_log, engine_log = log_dir / f"{stem}.gap.log", log_dir / f"{stem}.engine.log"
        records = []
        for task in batch_tasks:
            key = family_key(task)
            records.append(gap_record(task, tested_indices.get(key, []), safe_indices.get(key, [])))
        task_path.write_text("ASX_TASKS := [\n" + ",\n".join(records) + "\n];\n", encoding="utf-8")

        with input_path.open("w", encoding="utf-8") as output, gap_log.open("w", encoding="utf-8") as diagnostics:
            gap = subprocess.run(
                ["gap", "-q", "-b", "-c", f'Read("{task_path}"); Read("{run_exporter}");'],
                cwd=ROOT, stdin=subprocess.DEVNULL, text=True,
                stdout=output, stderr=diagnostics,
            )
        if gap.returncode != 0:
            raise RuntimeError(f"GAP failed for {stem}; see {gap_log}")
        exported_lines = [line.strip() for line in input_path.read_text().splitlines() if line.strip()]
        header = exported_lines[0] if exported_lines else ""
        exported = sum(line == "action" for line in exported_lines[1:])
        if header != "AFFINE_SAXL_V1":
            raise RuntimeError(f"bad GAP export header in {input_path}: {header!r}")

        rows: list[dict[str, object]] = []
        engine_jobs_used = 0
        if exported:
            engine_jobs_used = run_engine_parallel(
                engine, input_path, result_path, engine_log,
                args.engine_jobs, args.pair_budget,
            )
            rows = [json.loads(line) for line in result_path.read_text().splitlines() if line.strip()]
            if len(rows) != exported:
                raise RuntimeError(f"{stem}: exported {exported} actions but got {len(rows)} results")
        else:
            result_path.write_text("")
            engine_log.write_text("no uncovered catalogue entries\n")
            state["empty_rounds"] += len(batch_tasks)

        unsafe_families: set[str] = set()
        batch_labels: set[str] = set()
        for row in rows:
            label = str(row["label"])
            if label in seen_labels or label in batch_labels:
                raise RuntimeError(f"round exporter repeated tested label {label}")
            batch_labels.add(label)
            status = str(row["status"])
            safe = status in {
                "density_certificate", "fixed_space_union_bound", "two_sum_certificate"
            } or (row.get("base_size") == "2" and row.get("diameter") in (1, 2))
            expected = safe or (
                row.get("base_size") == ">2" and status == "no_regular_vector"
            ) or status == "diameter_at_least_3"
            if status in {"pair_budget_exceeded", "diameter_cap_reached"}:
                raise RuntimeError(f"resource guard hit for {label}: {row}")
            if not expected:
                raise RuntimeError(f"unexpected result for {label}: {row}")

        with master_results.open("a") as master:
            for row in rows:
                label = str(row["label"])
                seen_labels.add(label)
                task_data, k = task_from_label(label)
                key = family_key(task_data)
                tested_indices.setdefault(key, []).append(k)
                state["nodes"] += 1
                master.write(json.dumps(row, sort_keys=True) + "\n")
                status = str(row["status"])
                if status in {"density_certificate", "fixed_space_union_bound",
                              "two_sum_certificate"} or (
                    row.get("base_size") == "2" and row.get("diameter") in (1, 2)
                ):
                    state["safe_nodes"] += 1
                    safe_indices.setdefault(key, []).append(k)
                elif row.get("base_size") == ">2" and status == "no_regular_vector":
                    state["no_base_two_nodes"] += 1
                    unsafe_families.add(key)
                elif status == "diameter_at_least_3":
                    state["candidate_nodes"] += 1
                    unsafe_families.add(key)
                    with candidates_path.open("a") as candidate_file:
                        candidate_file.write(json.dumps(row, sort_keys=True) + "\n")
                else:
                    raise AssertionError(f"prevalidated result became unexpected: {row}")

        for key in tested_indices:
            tested_indices[key] = sorted(set(tested_indices[key]))
        for key in safe_indices:
            safe_indices[key] = sorted(set(safe_indices[key]))
        for task in batch_tasks:
            if family_key(task) in unsafe_families:
                next_task = dict(task)
                next_task["round"] = int(task["round"]) + 1
                state["queue"].append(next_task)

        state["processed_task_keys"] = sorted(processed)
        state["seen_node_labels"] = sorted(seen_labels)
        state["tested_indices"], state["safe_indices"] = tested_indices, safe_indices
        state["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        atomic_json(state_path, state)
        shard = {
            "schema": "IRREDSOL_CATALOGUE_ROUND_SHARD_V1", "batch": state["batch"],
            "tasks": len(batch_tasks), "exported_nodes": exported,
            "families_advanced": len(unsafe_families),
            "engine_jobs": engine_jobs_used,
            "task_sha256": sha256(task_path), "input_sha256": sha256(input_path),
            "result_sha256": sha256(result_path),
        }
        with shard_manifest.open("a") as destination:
            destination.write(json.dumps(shard, sort_keys=True) + "\n")
        print(
            f"{stem}: tasks={len(batch_tasks)} nodes={exported} "
            f"advanced={len(unsafe_families)} queue={len(state['queue'])} "
            f"safe={state['safe_nodes']} no_base2={state['no_base_two_nodes']} "
            f"candidates={state['candidate_nodes']}", flush=True,
        )

    state["processed_task_keys"] = sorted(processed)
    state["seen_node_labels"] = sorted(seen_labels)
    state["tested_indices"], state["safe_indices"] = tested_indices, safe_indices
    if not state["queue"]:
        state["completed_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        state["status"] = "complete"
    else:
        state["status"] = "paused"
    atomic_json(state_path, state)
    print(
        f"CATALOGUE_ROUNDS_{state['status'].upper()} batches={state['batch']} "
        f"nodes={state['nodes']} queue={len(state['queue'])} "
        f"candidates={state['candidate_nodes']} run_dir={run_dir}", flush=True,
    )


if __name__ == "__main__":
    main()
