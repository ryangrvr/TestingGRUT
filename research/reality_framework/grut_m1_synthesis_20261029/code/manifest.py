"""Create or verify SHA-256 hashes for the packet deliverables."""
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def payload():
 files={}
 for p in sorted(ROOT.rglob('*')):
  if not p.is_file() or p.name=='MANIFEST.json' or any((part in ('.work','__pycache__') or part.startswith('.work_')) for part in p.relative_to(ROOT).parts):continue
  files[str(p.relative_to(ROOT))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 return {'algorithm':'sha256','files':files}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();out=payload();path=ROOT/'MANIFEST.json'
if args.check:
 if out!=json.loads(path.read_text()):raise SystemExit('FAIL packet manifest mismatch')
 print('PASS manifest:',len(out['files']),'files')
else:
 path.write_text(json.dumps(out,indent=2)+'\n');print('Manifest written:',len(out['files']),'files')
