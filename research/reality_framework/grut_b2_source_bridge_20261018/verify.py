#!/usr/bin/env python3
"""Verify archived hashes and immutable primary/replay control agreement."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = json.loads((root / "MANIFEST.json").read_text())
for name, expected in manifest["files"].items():
    assert digest(root / name) == expected, name
pins = json.loads((root / "EXECUTED_SOURCES.json").read_text())
for name, expected in pins["files"].items():
    assert digest(root / name) == expected, "executed source " + name
primary = json.loads((root / "PRIMARY_WORKER.json").read_text())
replay = json.loads((root / "REPLAY.json").read_text())
actual = digest(root / "CONTROLS.json")
assert primary["output_sha256"] == replay["output_sha256"] == actual
assert primary["exit_code"] == replay["exit_code"] == 0
assert primary["wall_seconds"] + replay["wall_seconds"] < 60
assert max(primary["max_rss_kib"], replay["max_rss_kib"]) < 1024**2
assert sum(p.stat().st_size for p in root.iterdir() if p.is_file()) < 25*1024**2
results = json.loads((root / "CONTROLS.json").read_text())
assert results["status"] == "all_frozen_controls_pass"
assert results["witness"]["exact_return_modulus"] == 1
print(json.dumps({"status": "verified", "manifest_files": len(manifest["files"]),
                  "checks": len(results["checks"]), "output_sha256": actual}))
