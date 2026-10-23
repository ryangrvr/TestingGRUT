"""Pin packet files and consumed parent inputs; validate release integrity."""
import argparse,hashlib,json,py_compile,re
from pathlib import Path
import numpy as np
from model import ROOT,PARENT


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def validate():
 for path in ROOT.rglob('*.py'):py_compile.compile(str(path),doraise=True)
 def no_constant(value):raise ValueError('Non-finite JSON constant: '+value)
 for path in ROOT.glob('*.json'):json.loads(path.read_text(),parse_constant=no_constant)
 for path in (ROOT/'data').glob('*.csv'):
  assert np.isfinite(np.loadtxt(path,delimiter=',',skiprows=1)).all(),path
 for path in (ROOT/'data').glob('*.npz'):
  with np.load(path) as a:
   for key in a.files:assert np.isfinite(a[key]).all(),(path,key)
 for path in ROOT.glob('*.md'):
  for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
   if '://' not in target and not target.startswith('#'):
    assert (path.parent/target.split('#')[0]).exists(),(path.name,target)


def run(check=False):
 validate()
 files=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json']
 consumed=[PARENT/'RESULTS.json',*[PARENT/'data'/f'profile_Q{q}.csv' for q in [19972,39944,79988]]]
 current={'parent_commit':'fa3bc97b52f44bdd625da8b40fdc7015837517a7',
  'files':{str(p.relative_to(ROOT)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in files},
  'consumed_parent_files':{str(p.relative_to(ROOT.parent)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in consumed}}
 manifest=ROOT/'MANIFEST.json'
 if check:assert json.loads(manifest.read_text())==current,'Manifest differs from present packet or inputs'
 else:manifest.write_text(json.dumps(current,indent=2)+'\n')
 print(json.dumps({'manifest':'verified' if check else 'written','packet_files':len(files),'parent_inputs':len(consumed),'syntax_finite_arrays_and_local_links':'pass'}))


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');run(p.parse_args().check)
