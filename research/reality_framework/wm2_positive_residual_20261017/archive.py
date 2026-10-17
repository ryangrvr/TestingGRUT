"""Seal local deliverables and immutable imported-source dependencies."""
from pathlib import Path
import hashlib
import json

BASE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=sorted([p for p in BASE.iterdir() if p.is_file() and p.name!='MANIFEST.json']+
        list((BASE/'results').glob('*.json')))
    deps=['wm2_history_formation_20261016/'+x for x in ['active.py','controls.py','manifest.json']]
    deps+=['wm2_prepared_disorder_20261016/model.py','wm2_conditioned_receiver_20261016/manifest.json']
    result=dict(parent='b4c79020e828b7876d54bd68e120237fc147396e',
        charter_commit='4bf14483a0804980244765906e939e4dbaf29834',
        sha256={str(p.relative_to(BASE)):sha(p) for p in files},
        dependencies={name:sha(BASE.parent/name) for name in deps})
    (BASE/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
