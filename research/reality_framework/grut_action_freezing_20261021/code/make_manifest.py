from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parents[1]
files={str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(p.rglob('*')) if f.is_file() and f.name!='MANIFEST.json'}
(p/'MANIFEST.json').write_text(json.dumps(dict(algorithm='sha256',files=files),indent=2,sort_keys=True)+'\n')
print('Manifest:',len(files),'exact payload files')
