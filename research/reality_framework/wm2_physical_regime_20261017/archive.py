"""Hash owned artifacts and every frozen runtime/control dependency."""
from pathlib import Path
import hashlib
import json

BASE=Path(__file__).resolve().parent
S8=BASE.parent/'wm2_schmidt_receiver_20261017'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    files=sorted([p for p in BASE.iterdir() if p.is_file() and p.name!='MANIFEST.json' and p.suffix!='.partial']+
        list((BASE/'results').glob('*.json')))
    names=list(json.loads((S8/'MANIFEST.json').read_text())['dependencies'])
    names+=['wm2_schmidt_receiver_20261017/'+x for x in ['MANIFEST.json','bounds.py','gates.py',
        'allocation.py','controls.py','sweep_control.py','mps.py','mps_global.py','environment.json']]
    result=dict(parent='cb361f8f230371515f1806ab50ca62657a3088be',
        freeze='85fc1bd2c4592a57c814bc5b2879aab65fdb85aa',
        sha256={str(p.relative_to(BASE)):sha(p) for p in files},
        dependencies={name:sha(BASE.parent/name) for name in names})
    (BASE/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':run()
