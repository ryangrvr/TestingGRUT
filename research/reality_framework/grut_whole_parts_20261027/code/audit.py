"""Validate this small packet and its source pins, without replaying old physics."""
import ast
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from common import ROOT, save


def main():
    sources = json.loads((ROOT/"SOURCES.json").read_text())
    repo = ROOT.parents[2]
    for row in sources["local_inputs_and_context"]:
        actual = hashlib.sha256((repo/row["path"]).read_bytes()).hexdigest()
        assert actual == row["sha256"], row["path"]
    manifest = json.loads((ROOT/"MANIFEST.json").read_text())
    for row in manifest["files"]:
        p = ROOT/row["path"]
        assert p.stat().st_size == row["bytes"], str(p)
        assert hashlib.sha256(p.read_bytes()).hexdigest() == row["sha256"], str(p)
    known = {row["path"] for row in manifest["files"]} | {"MANIFEST.json", "AUDIT.json"}
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*")
              if p.is_file() and "__pycache__" not in p.parts}
    assert actual == known, {"missing": list(known-actual), "extra": list(actual-known)}
    json_count = python_count = link_count = 0
    for p in ROOT.rglob("*"):
        if not p.is_file() or "__pycache__" in p.parts:
            continue
        if p.suffix == ".json":
            obj = json.loads(p.read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
            json_count += 1
        elif p.suffix == ".py":
            ast.parse(p.read_text())
            python_count += 1
        elif p.suffix == ".md":
            for target in re.findall(r'\]\(([^)]+)\)', p.read_text()):
                if "://" in target or target.startswith("#"):
                    continue
                assert (p.parent/target.split("#")[0]).exists(), (str(p), target)
                link_count += 1
    with np.load(ROOT/"HISTORY_MODEL.npz") as a:
        for key in a.files:
            assert np.isfinite(a[key]).all(), key
    image = (ROOT/"WHOLE_PARTS.png").read_bytes()
    assert image[:8] == b'\x89PNG\r\n\x1a\n'
    verification = json.loads((ROOT/"VERIFICATION.json").read_text())
    assert verification["all_pass"] and len(verification["groups"]) == 8
    save("AUDIT.json", {"all_pass": True, "source_pins": len(sources["local_inputs_and_context"]),
                        "manifest_hashed_files": len(manifest["files"]), "json_files": json_count,
                        "python_files": python_count, "local_links": link_count,
                        "finite_model_arrays": True, "mathematical_groups": 8,
                        "independent_review": False})
    print(json.dumps({"audit": "pass", "source_pins": len(sources["local_inputs_and_context"]),
                      "mathematical_groups": 8}))


if __name__ == "__main__":
    main()
