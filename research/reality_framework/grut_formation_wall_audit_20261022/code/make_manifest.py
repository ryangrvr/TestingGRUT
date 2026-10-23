"""Integrity/syntax checks and exact hashes, separate from physical verification."""
import argparse,hashlib,json,py_compile,re
import numpy as np
from calculate import ROOT,EMP


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run(check=False):
 for p in ROOT.rglob('*.py'):py_compile.compile(str(p),doraise=True)
 def nonfinite(value):raise ValueError(value)
 for p in ROOT.glob('*.json'):json.loads(p.read_text(),parse_constant=nonfinite)
 for p in (ROOT/'data').glob('*.csv'):assert np.isfinite(np.loadtxt(p,delimiter=',',skiprows=1)).all(),p
 for p in ROOT.glob('*.md'):
  prose=re.sub(r'\\\[.*?\\\]','',p.read_text(),flags=re.S)
  for target in re.findall(r'\]\(([^)\n]+)\)',prose):
   if '://' not in target and not target.startswith('#'):assert (p.parent/target.split('#')[0]).exists(),(p,target)
 files=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json']
 inputs=[EMP/'DETECTOR_RESULTS.json',EMP/'HISTORY_RESULTS.json']
 content={'parent_commit':'8b0d96b219ea6616de23d0adbb67d4dc8be511ea',
  'files':{str(p.relative_to(ROOT)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in files},
  'consumed_parent_files':{str(p.relative_to(ROOT.parent)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in inputs}}
 path=ROOT/'MANIFEST.json'
 if check:assert json.loads(path.read_text())==content
 else:path.write_text(json.dumps(content,indent=2)+'\n')
 print(json.dumps({'manifest':'verified' if check else 'written','packet_files':len(files),'parent_inputs':len(inputs),'syntax_finite_values_local_links':'pass'}))


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');run(p.parse_args().check)
