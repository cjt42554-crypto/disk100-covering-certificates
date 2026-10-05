"""Portable exact centerwise transport; this new adapter is unexecuted."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "common"))
from exact_radial import centers_at, read_npy
from public_inputs import load_manifest, read_bound


def main():
    manifest = load_manifest()
    for source in ["transport/source/verify_transport.py", "common/exact_radial.py", "common/public_inputs.py"]:
        read_bound(source, manifest)
    mapping = json.loads(read_bound("transport/INPUTS.json", manifest))
    def load(item):
        return read_bound(item["file"], manifest, item["sha256"])
    arrays = []
    for key in ["source", "ordinary_best"]:
        payload, endian = read_npy(load(mapping[key]), (100, 2))
        arrays.append(centers_at(payload, endian))
    inherited = json.loads(load(mapping["inherited_transport_result"]))
    delta = F(1, 100000000)
    for i, (a, b) in enumerate(zip(*arrays)):
        differences = [F.from_float(b[j]) - F.from_float(a[j]) for j in range(2)]
        distance_squared = sum(x * x for x in differences)
        if distance_squared > delta * delta:
            raise ValueError("centerwise displacement exceeds delta")
        saved = inherited["rows"][i]
        if saved["center"] != i:
            raise ValueError("center index mismatch")
        if saved["difference_xy"] != [[x.numerator, x.denominator] for x in differences]:
            raise ValueError("saved exact coordinate difference mismatch")
        if saved["distance_squared"] != [distance_squared.numerator, distance_squared.denominator]:
            raise ValueError("saved squared-distance mismatch")
        if saved["distance_squared_le_delta_squared"] is not True:
            raise ValueError("saved transport acceptance mismatch")
    source_cover = F(14550606453, 125000000000)
    transported_cover = source_cover + delta
    source_exclusion = F(2362504, 1000000000)
    transported_exclusion = source_exclusion - delta
    if transported_cover != F(*inherited["exact_closed_D_cover_radius_fraction"]):
        raise ValueError("transport radius mismatch")
    print(json.dumps({"status": "PASS_EXACT_CENTERWISE_DISPLACEMENT",
                      "all100_checks": 100, "delta": [1, 100000000],
                      "conditional_transport_cover_radius": [transported_cover.numerator, transported_cover.denominator],
                      "conditional_finite_reference_strict_distance_bound": [transported_exclusion.numerator, transported_exclusion.denominator],
                      "prerequisites": "The source exact disk cover and source finite Hall certificate must also hold. This verifier checks displacement, not their separate proofs."}, indent=2))


if __name__ == "__main__":
    main()
