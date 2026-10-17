from pathlib import Path
import hashlib
import json

BASE=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=sorted([p for p in BASE.iterdir() if p.is_file() and p.name!='MANIFEST.json' and p.suffix!='.partial']+
        list((BASE/'results').glob('*.json')))
    deps=['wm2_positive_residual_20261017/'+x for x in ['native.cpp','native.py','enclosures.py','rotation.py','layout.py','solver.py','worker.py','MANIFEST.json']]
    deps+=['wm2_history_formation_20261016/active.py','wm2_history_formation_20261016/manifest.json',
        'wm2_prepared_disorder_20261016/model.py']
    result=dict(parent='a1dd220ee28cf8fd985a7c144922e0bb6c299902',
        original_charter='0e42a5d3023c7ec45029f38ef918c47f62b6c603',
        G1_charter='b15fcf0618ac53e2a5a69f216ac346e16fa3f6ff',
        sha256={str(p.relative_to(BASE)):sha(p) for p in files},
        dependencies={name:sha(BASE.parent/name) for name in deps})
    (BASE/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
