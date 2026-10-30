"""Verify exact inherited inputs and package the additive research packet."""
from pathlib import Path
import json,hashlib,re,ast
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[2]
def main():
 pins=json.loads((ROOT/'INPUT_PINS.json').read_text());errors=[]
 for name,item in pins['files'].items():
  data=(REPO/name).read_bytes();blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
  if blob!=item['git_blob_sha1'] or hashlib.sha256(data).hexdigest()!=item['sha256']:errors.append('Input mismatch: '+name)
 files=[x for x in ROOT.rglob('*') if x.is_file() and not any(t in x.relative_to(ROOT).parts for t in ['.work','__pycache__']) and x.name!='MANIFEST.json'];links=0
 for path in files:
  if path.suffix=='.py':ast.parse(path.read_text())
  if path.suffix=='.json':json.loads(path.read_text())
  if path.suffix=='.npz':
   with np.load(path) as data:
    if not all(np.isfinite(data[k]).all() for k in data.files):errors.append('Nonfinite array: '+path.name)
  if path.suffix=='.md':
   for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
    if re.match(r'\w+://',target) or target.startswith('#'):continue
    links+=1
    if not (path.parent/target.split('#')[0]).exists():errors.append('Broken link: '+str(path)+': '+target)
 v=json.loads((ROOT/'VERIFICATION.json').read_text());n=json.loads((ROOT/'FOLLOWUP_NOTES_CHECK.json').read_text())
 if not v['all_pass'] or n['passed']!=n['total']:errors.append('Unpassed computational check')
 manifest={'parent_commit':pins['parent_commit'],'parent_tree':pins['parent_tree'],'input_pins_checked':len(pins['files']),'local_links_checked':links,'same_author_dynamics_checks':v['passed'],'followup_notes_checks':n['passed'],'files':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)},'scope':'Manifest of additive artifacts; excludes itself, .work, and Python caches. Check counts are implementation/counterexample controls, not external review.'}
 if errors:raise RuntimeError('\n'.join(errors))
 (ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'files':len(files),'input_pins':len(pins['files']),'local_links':links,'checks':v['passed']+n['passed'],'errors':errors}))
if __name__=='__main__':main()
