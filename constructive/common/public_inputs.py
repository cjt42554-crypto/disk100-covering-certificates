"""Load project-relative mathematical files bound by the clean manifest."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_manifest():
    return json.loads((ROOT / "MANIFEST_SHA256.json").read_text())


def read_bound(relative_path, manifest, expected_sha256=None):
    path = (ROOT / relative_path).resolve()
    if not path.is_relative_to(ROOT.resolve()):
        raise ValueError("input leaves constructive project")
    canonical = path.relative_to(ROOT.resolve()).as_posix()
    if canonical != relative_path:
        raise ValueError("input path is not canonical project-relative form")
    data = path.read_bytes()
    actual = hashlib.sha256(data).hexdigest()
    if manifest[relative_path] != actual:
        raise ValueError("manifest hash mismatch: " + relative_path)
    if expected_sha256 is not None and expected_sha256 != actual:
        raise ValueError("input mapping hash mismatch: " + relative_path)
    return data
