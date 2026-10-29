"""Seal this packet's non-cache files and check pinned inherited inputs."""
from setup import *
import ast,re

def main():
    pins=json.loads((ROOT/'INPUT_PINS.json').read_text())['sha256']
    changed=[name for name,digest in pins.items() if sha(RF/name)!=digest]
    if changed:raise RuntimeError('Inherited inputs changed: '+repr(changed))
    files=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '.work' in p.parts or '__pycache__' in p.parts:continue
        if p.name in ('MANIFEST.json','AUDIT.json'):continue
        files.append(p)
        if p.suffix=='.py':ast.parse(p.read_text())
        if p.suffix=='.json':json.loads(p.read_text())
    broken=[]
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if not link.startswith(('http:','https:','#')) and not (p.parent/link.split('#')[0]).exists():broken.append((p.name,link))
    if broken:raise RuntimeError('Broken local links: '+repr(broken))
    write_json('MANIFEST.json',{'sha256':{str(p.relative_to(ROOT)):sha(p) for p in files}})
    write_json('AUDIT.json',{'input_pins_unchanged':True,'pinned_inputs':len(pins),'sealed_files':len(files),
        'local_links_valid':True,'finite_arrays_checked_in':'VERIFICATION.json','independent_author_review':False})
    print('Sealed',len(files),'files; pins and local links pass.')

if __name__=='__main__':main()
