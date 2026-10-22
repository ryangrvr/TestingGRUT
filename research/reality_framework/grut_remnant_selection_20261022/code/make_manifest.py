#!/usr/bin/env python3
"""Finalize hashes; verifier never calls this script."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
       for p in sorted(ROOT.rglob('*'))
       if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
(ROOT/'MANIFEST.json').write_text(json.dumps({'algorithm':'sha256','files':files},indent=2)+'\n')
print(f'Hashed {len(files)} files; MANIFEST.json excludes itself.')
