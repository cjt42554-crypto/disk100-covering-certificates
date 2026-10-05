"""One fresh exact finite reference certificate, launched only by root.

Floating guidance selects a rational threshold. It never accepts a witness.
No geometric cover score, optimizer, pose fit or ordinary restore is called.
"""
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import struct
import sys
import time

from exact_radial import (centers_at, coordinate_digest, exact_norm_squares,
                         hall_accept, read_npy, sha256)

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT.parent
DENOMINATOR = 1_000_000_000
EPSILON_NUMERATOR = 2000  # 2000/10^9 = 1/500000 exactly.


def write_new(path, value):
    with path.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def utc():
    return datetime.now(timezone.utc).isoformat()


def admitted_inputs():
    plan_bytes = (ROOT / "PLAN.json").read_bytes()
    plan = json.loads(plan_bytes)
    launch = json.loads((ROOT / "EXECUTION_STARTED.json").read_text())
    if launch["plan_sha256"] != sha256(plan_bytes):
        raise ValueError("plan changed after one-shot launch")
    if plan["status"] != "READY_FOR_ONE_ROOT_EXECUTION":
        raise ValueError("plan not ready")
    if plan["epsilon"] != {"numerator": 1, "denominator": 500000}:
        raise ValueError("epsilon changed")
    if plan["reference_count"] != 35599 or plan["centers_per_array"] != 100:
        raise ValueError("finite/all100 scope changed")
    for name, expected in plan["source_hashes"].items():
        if sha256((ROOT / name).read_bytes()) != expected:
            raise ValueError("source changed: " + name)
    binding_bytes = (ROOT / "INPUT_BINDING.json").read_bytes()
    if sha256(binding_bytes) != plan["input_binding_sha256"]:
        raise ValueError("input binding changed")
    binding = json.loads(binding_bytes)
    if binding["status"] != "ALL_PINNED_BYTES_AND_NPY_HEADERS_BOUND_NO_PROOF_EXECUTION":
        raise ValueError("input binding not ready")
    # Every relevant input is read once, hash-checked, then reused from memory.
    payloads = {}
    for name, item in binding["inputs"].items():
        data = (BASE / name).read_bytes()
        if sha256(data) != item["sha256"] or len(data) != item["bytes"]:
            raise ValueError("pinned input changed: " + name)
        payloads[name] = data
    manifest = json.loads(payloads[binding["reference_manifest"]])
    inventory = json.loads(payloads[binding["reference_inventory"]])
    reference_prefix = "radial_guarded_endpoint_audit/"
    for name, item in binding["inputs"].items():
        if name.startswith(reference_prefix) and name != binding["reference_manifest"]:
            if manifest[name[len(reference_prefix):]] != item["sha256"]:
                raise ValueError("reference manifest does not bind " + name)
    if inventory["unique_unordered100_arrays"] != 35599:
        raise ValueError("reference count mismatch")
    if len(inventory["entries"]) != 35599 or len(inventory["chunks"]) != 36:
        raise ValueError("entry/chunk count mismatch")
    if payloads[binding["current_source"]] != payloads[binding["sealed_endpoint"]]:
        raise ValueError("actual source differs from sealed endpoint")
    return plan, binding, inventory, payloads


def main():
    if not __debug__ or not sys.flags.no_site or not sys.flags.dont_write_bytecode:
        raise RuntimeError("root must launch python3 -S -B")
    if sys.argv[1:] != ["--from-once-wrapper"]:
        raise RuntimeError("use source/execute_once.py")
    plan, binding, inventory, payloads = admitted_inputs()
    results_dir = ROOT / "results"
    results_dir.mkdir()  # refuses all existing result directories, even on failure.
    begin, tick = utc(), time.monotonic()
    processed = accepted = buffered_accepted = 0
    minimum_buffer = None
    minimum_guidance_gap = None
    current_payload, current_endian = read_npy(
        payloads[binding["current_source"]], (100, 2))
    current = centers_at(current_payload, current_endian)
    current_exact = exact_norm_squares(current)
    current_norms = sorted(math.hypot(x, y) for x, y in current)
    norm_payload, norm_endian = read_npy(
        payloads[binding["norm_index"]], (35599, 100))
    chunk_payloads = {}
    for chunk in inventory["chunks"]:
        name = "radial_guarded_endpoint_audit/" + chunk["file"]
        if binding["inputs"][name]["sha256"] != chunk["file_sha256"]:
            raise ValueError("inventory/chunk hash mismatch")
        chunk_payloads[chunk["file"]] = read_npy(
            payloads[name], (chunk["arrays"], 100, 2))
    visited = set()
    witness_path = results_dir / "WITNESSES.jsonl"
    with witness_path.open("x", encoding="utf-8") as witnesses:
        for index, entry in enumerate(inventory["entries"]):
            chunk_name, row_index = entry["chunk"], entry["row"]
            pair = (chunk_name, row_index)
            if pair in visited:
                raise ValueError("reference row repeated")
            visited.add(pair)
            chunk, endian = chunk_payloads[chunk_name]
            reference = centers_at(chunk, endian, row_index)
            digest = coordinate_digest(reference)
            if digest != entry["unordered_coordinate_sha256"]:
                raise ValueError("reference actual-coordinate digest mismatch")
            reference_exact = exact_norm_squares(reference)
            guide = struct.unpack_from(norm_endian + "100d", norm_payload, index * 800)
            if not all(math.isfinite(x) and x >= 0 for x in guide):
                raise ValueError("nonfinite/negative numeric guidance")
            if any(guide[i] > guide[i + 1] for i in range(99)):
                raise ValueError("numeric guidance is not sorted")
            # This selector may be arbitrarily inaccurate: exact acceptance below
            # is what establishes the theorem. No sorting of exact roots is used.
            k = max(range(100), key=lambda i: abs(current_norms[i] - guide[i]))
            gap = abs(current_norms[k] - guide[k])
            minimum_guidance_gap = gap if minimum_guidance_gap is None else min(minimum_guidance_gap, gap)
            low_current = current_norms[k] < guide[k]
            lo, hi = sorted((current_norms[k], guide[k]))
            epsilon_float = EPSILON_NUMERATOR / DENOMINATOR
            threshold = math.floor(((lo + hi - epsilon_float) / 2) * DENOMINATOR)
            buffer = max(0, math.floor((gap - epsilon_float) * DENOMINATOR / 4) - 4)
            low, high = (current_exact, reference_exact) if low_current else (reference_exact, current_exact)
            witness = hall_accept(low, high, threshold, DENOMINATOR, EPSILON_NUMERATOR)
            # Optional exact radial margin: a second Hall pair at t-buffer and
            # t+epsilon+buffer proves the stronger radius epsilon+2*buffer.
            buffered = hall_accept(low, high, threshold - buffer,
                                   DENOMINATOR, EPSILON_NUMERATOR + 2 * buffer)
            if not witness["accepted"] or not buffered["accepted"]:
                raise ValueError("integer Hall deficit absent for reference " + str(index))
            accepted += 1
            buffered_accepted += 1
            minimum_buffer = buffer if minimum_buffer is None else min(minimum_buffer, buffer)
            record = {
                "reference_index": index, "chunk": chunk_name, "row": row_index,
                "reference_unordered_coordinate_sha256": digest,
                "low_side": "current" if low_current else "reference",
                "threshold_numerator": threshold, "threshold_denominator": DENOMINATOR,
                "epsilon_numerator": EPSILON_NUMERATOR,
                "guidance_sorted_rank_0based": k,
                "guidance_gap_numeric_only": gap,
                "hall": witness, "buffer_numerator": buffer,
                "buffered_hall": buffered,
            }
            witnesses.write(json.dumps(record, separators=(",", ":"), allow_nan=False) + "\n")
            processed += 1
            if processed % 1000 == 0 or processed == 35599:
                print(json.dumps({"stage": "EXACT_INTEGER_RADIAL_HALL", "references": processed,
                                  "accepted": accepted, "all100": True}), flush=True)
    expected_rows = {(c["file"], i) for c in inventory["chunks"] for i in range(c["arrays"])}
    if visited != expected_rows or len(visited) != 35599:
        raise ValueError("reference chunks not covered exactly once")
    result = {
        "status": "PASS_EXACT_FINITE_ORIGIN_FIXED_O2_ALL100_EXCLUSION",
        "reference_count": processed, "hall_witnesses_accepted": accepted,
        "buffered_hall_witnesses_accepted": buffered_accepted,
        "centers_per_array": 100,
        "epsilon": {"numerator": 1, "denominator": 500000},
        "conclusion": "For each pinned finite reference, every origin-fixed orthogonal transform and every all100 bijection has maximum distance strictly greater than epsilon.",
        "current_source_sha256": binding["inputs"][binding["current_source"]]["sha256"],
        "current_ordered_coordinate_sha256": coordinate_digest(current),
        "reference_manifest_sha256": binding["inputs"][binding["reference_manifest"]]["sha256"],
        "reference_inventory_sha256": binding["inputs"][binding["reference_inventory"]]["sha256"],
        "witnesses_file": "results/WITNESSES.jsonl", "witnesses_sha256": sha256(witness_path.read_bytes()),
        "minimum_exact_buffer": {"numerator": minimum_buffer, "denominator": DENOMINATOR},
        "stronger_common_certified_epsilon": {"numerator": EPSILON_NUMERATOR + 2 * minimum_buffer, "denominator": DENOMINATOR},
        "minimum_selected_norm_gap_numeric_only": minimum_guidance_gap,
        "integer_coordinate_conversions": 2 * 100 * (35599 + 1),
        "integer_count_masks": 4 * 35599,
        "numeric_guidance_only_no_numeric_acceptance": True,
        "translation_allowed": False, "scaling_allowed": False,
        "references_are_finite_saved_arrays_only": True,
        "exhaustive_historical_inventory_claim": False,
        "basin_exclusion_claim": False, "global_optimum_claim": False,
        "new_cover_score_calls": 0, "optimizer_calls": 0, "original_calls": 0,
        "plan_sha256": sha256((ROOT / "PLAN.json").read_bytes()),
        "utc_start": begin, "utc_end": utc(), "elapsed_seconds": time.monotonic() - tick,
    }
    write_new(ROOT / "RESULT.json", result)
    print(json.dumps({"status": result["status"], "references": processed,
                      "elapsed_seconds": result["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
