#!/usr/bin/env python3
"""finalize.py -- write MANIFEST.json (sha256 of every packet file except the manifest itself)."""
import hashlib
import json
import sys
from pathlib import Path

PKT = Path(__file__).resolve().parent.parent
files = {}
for p in sorted(PKT.rglob("*")):
    if p.is_file() and p.name != "MANIFEST.json" and "__pycache__" not in p.parts:
        files[str(p.relative_to(PKT))] = hashlib.sha256(p.read_bytes()).hexdigest()
man = {"packet": PKT.name, "generated_by": "code/finalize.py", "n_files": len(files), "files": files,
       "python": sys.version.split()[0]}
(PKT / "MANIFEST.json").write_text(json.dumps(man, indent=1) + "\n")
print(f"MANIFEST.json written: {len(files)} files")
