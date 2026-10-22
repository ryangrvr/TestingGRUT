"""Refresh only this proposal's content manifest."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
files = {str(p.relative_to(root)): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size}
         for p in sorted(root.rglob("*")) if p.is_file() and p.name != "MANIFEST.json" and "__pycache__" not in p.parts}
(root/"MANIFEST.json").write_text(json.dumps({"files": files}, indent=2)+"\n")
print(len(files), "files plus manifest")
