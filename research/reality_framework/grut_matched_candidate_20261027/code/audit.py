"""Validate files and pinned imports without rerunning numerical production."""
import ast,hashlib,json,re
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[2]
s=json.loads((ROOT/'SOURCES.json').read_text())
for row in s['local_sources']:
    assert hashlib.sha256((REPO/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
m=json.loads((ROOT/'MANIFEST.json').read_text())
for row in m['files']:
    b=(ROOT/row['path']).read_bytes()
    assert len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],row['path']
known={row['path'] for row in m['files']}|{'MANIFEST.json','AUDIT.json'}
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
assert actual==known
links=0
for p in ROOT.rglob('*'):
    if not p.is_file() or '__pycache__' in p.parts:continue
    if p.suffix=='.py':ast.parse(p.read_text())
    elif p.suffix=='.json':json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
    elif p.suffix=='.csv':assert np.isfinite(np.loadtxt(p,delimiter=',',skiprows=1)).all()
    elif p.suffix=='.npz':
        with np.load(p) as arr:assert all(np.isfinite(arr[k]).all() for k in arr.files)
    elif p.suffix=='.md':
        body=re.sub(r'\\\[.*?\\\]','',p.read_text(),flags=re.S)
        for link in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',body):
            if '://' in link or link.startswith('#'):continue
            assert (p.parent/link.split('#')[0]).exists(),link
            links+=1
v=json.loads((ROOT/'VERIFICATION.json').read_text());assert v['all_pass'] and len(v['groups'])==8
assert (ROOT/'M1.png').read_bytes()[:8]==b'\x89PNG\r\n\x1a\n'
result={'all_pass':True,'source_pins':len(s['local_sources']),'mathematical_groups':8,
        'manifest_files':len(m['files']),'local_links':links,'independent_review':False}
(ROOT/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
