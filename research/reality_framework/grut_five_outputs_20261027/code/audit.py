"""Check scope-limited packet integrity and pinned numerical inputs."""
import ast,hashlib,json,re
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[2]
sources=json.loads((ROOT/'SOURCES.json').read_text())
for row in sources['local_sources']:
    assert hashlib.sha256((REPO/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
manifest=json.loads((ROOT/'MANIFEST.json').read_text())
for row in manifest['files']:
    data=(ROOT/row['path']).read_bytes()
    assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],row['path']
known={row['path'] for row in manifest['files']}|{'MANIFEST.json','AUDIT.json'}
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
assert actual==known
links=0
for p in ROOT.rglob('*'):
    if not p.is_file() or '__pycache__' in p.parts:continue
    if p.suffix=='.json':json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
    elif p.suffix=='.py':ast.parse(p.read_text())
    elif p.suffix=='.npz':
        with np.load(p) as arrays:
            assert all(np.isfinite(arrays[k]).all() for k in arrays.files)
    elif p.suffix=='.md':
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' in link or link.startswith('#'):continue
            assert (p.parent/link.split('#')[0]).exists(),link
            links+=1
verification=json.loads((ROOT/'VERIFICATION.json').read_text())
assert verification['all_pass'] and len(verification['groups'])==6
assert (ROOT/'FIVE_OUTPUTS.png').read_bytes()[:8]==b'\x89PNG\r\n\x1a\n'
result={'all_pass':True,'source_pins':len(sources['local_sources']),'mathematical_groups':6,
        'manifest_files':len(manifest['files']),'local_links':links,'independent_review':False}
(ROOT/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
