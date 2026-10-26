"""Final source/data hashes; no temporary numerical caches or binaries."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
paths=[p for p in ROOT.rglob('*') if p.is_file() and p.suffix in {'.md','.py','.c','.csv','.json','.png'} and p.name!='MANIFEST.json' and '__pycache__' not in p.parts]
rows=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(paths)]
(ROOT/'MANIFEST.json').write_text(json.dumps(dict(files=rows,scope='Final additive research packet. Temporary external caches and bytecode excluded.'),indent=2)+'\n')
print('MANIFEST',len(rows),'files',sum(r['bytes'] for r in rows),'bytes')
