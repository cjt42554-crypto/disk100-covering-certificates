"""Portable exact verification of the saved finite radial Hall witnesses.

This public adapter has NOT been mathematically executed in preparation.
The original acceptance kernel is unchanged. No floating norm selector is needed.
"""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "common"))
from exact_radial import (centers_at, coordinate_digest, exact_norm_squares,
                         hall_accept, read_npy)
from public_inputs import load_manifest, read_bound


def main():
    manifest = load_manifest()
    for source in ["hall/source/check_all.py", "common/exact_radial.py", "common/public_inputs.py"]:
        read_bound(source, manifest)
    mapping = json.loads(read_bound("hall/INPUTS.json", manifest))
    def load(item):
        return read_bound(item["file"], manifest, item["sha256"])
    source_item = mapping["current_source"]
    payload, endian = read_npy(load(source_item), (100, 2))
    current = centers_at(payload, endian)
    if coordinate_digest(current) != mapping["current_ordered_coordinate_sha256"]:
        raise ValueError("actual current coordinate digest mismatch")
    current_exact = exact_norm_squares(current)
    inventory = json.loads(load(mapping["inventory"]))
    entries, chunks = inventory["entries"], inventory["chunks"]
    if len(entries) != 35599 or len(chunks) != 36 or mapping["reference_count"] != 35599:
        raise ValueError("finite reference scope mismatch")
    arrays = {}
    expected_rows = set()
    for item in chunks:
        arrays[item["file"]] = read_npy(load(item), (item["arrays"], 100, 2))
        expected_rows.update((item["file"], i) for i in range(item["arrays"]))
    rows = load(mapping["witnesses"]).decode("utf-8").splitlines()
    if len(rows) != len(entries):
        raise ValueError("witness count mismatch")
    visited = set()
    minimum_buffer = None
    for index, (entry, line) in enumerate(zip(entries, rows)):
        witness = json.loads(line)
        chunk = "hall/" + witness["chunk"]
        row = witness["row"]
        if witness["reference_index"] != index or (chunk, row) != (entry["chunk"], entry["row"]):
            raise ValueError("witness/inventory order mismatch")
        if (chunk, row) in visited:
            raise ValueError("duplicate reference row")
        visited.add((chunk, row))
        payload, endian = arrays[chunk]
        reference = centers_at(payload, endian, row)
        digest = coordinate_digest(reference)
        if digest != entry["unordered_coordinate_sha256"] or digest != witness["reference_unordered_coordinate_sha256"]:
            raise ValueError("reference actual IEEE coordinate identity mismatch")
        reference_exact = exact_norm_squares(reference)
        if witness["low_side"] == "current":
            low, high = current_exact, reference_exact
        elif witness["low_side"] == "reference":
            low, high = reference_exact, current_exact
        else:
            raise ValueError("invalid Hall direction")
        t, d, e, b = (witness[k] for k in ["threshold_numerator", "threshold_denominator", "epsilon_numerator", "buffer_numerator"])
        if not all(type(v) is int for v in (t, d, e, b)) or b < 0 or t < b:
            raise ValueError("invalid exact threshold or buffer")
        if d != 1000000000 or e != 2000:
            raise ValueError("epsilon must remain exactly 1/500000")
        basic = hall_accept(low, high, t, d, e)
        buffered = hall_accept(low, high, t - b, d, e + 2 * b)
        if basic != witness["hall"] or buffered != witness["buffered_hall"]:
            raise ValueError("saved count/mask differs from all100 integer recomputation")
        if not basic["accepted"] or not buffered["accepted"]:
            raise ValueError("strict integer Hall deficit absent")
        minimum_buffer = b if minimum_buffer is None else min(minimum_buffer, b)
    if visited != expected_rows or len(visited) != 35599:
        raise ValueError("reference chunk rows not covered exactly once")
    if 2000 + 2 * minimum_buffer != 2362504:
        raise ValueError("stronger common threshold differs from inherited result")
    print(json.dumps({"status": "PASS", "reference_count": len(visited),
                      "all100_indexed_centers": True,
                      "common_strict_distance_lower_bound": [2362504, 1000000000],
                      "scope": "Pinned finite references, origin-fixed O(2), every all100 bijection; no basin or global optimum claim."}, indent=2))


if __name__ == "__main__":
    main()
