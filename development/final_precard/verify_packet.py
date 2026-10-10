"""Verify the public review envelope without acquiring private reviewer inputs."""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]


def verify():
    p=ROOT/'development/final_precard/REVIEW_MANIFEST.json'
    raw=p.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=p.with_suffix('.json.sha256').read_text().strip():
        raise ValueError('Review manifest hash lock mismatch')
    m=json.loads(raw)
    if m['parent_sha']!='b89a0e607e71450fce07d5cc55daae8224acad18' or m['external_approval'] or m['card_1_authorized']:
        raise ValueError('Review envelope changed its authority')
    subprocess.run(['git','merge-base','--is-ancestor',m['work_branch_sha'],'HEAD'],cwd=ROOT,check=True,capture_output=True)
    for name,expected in m['public_files'].items():
        p=ROOT/name; data=p.read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected['sha256']:
            raise ValueError('Reviewed public bytes changed: '+name)
        original=subprocess.check_output(['git','show',m['work_branch_sha']+':'+name],cwd=ROOT)
        if original!=data: raise ValueError('Source pin differs from manifest bytes: '+name)
    r=ROOT/'development/final_precard/evidence/RESULTS.json'
    if hashlib.sha256(r.read_bytes()).hexdigest()!=m['results_sha256']:
        raise ValueError('Results digest mismatch')
    result=json.loads(r.read_bytes())
    if result['status']!='READY FOR EXTERNAL CHECK' or result['remaining_semantic_mismatches'] or result['external_approval']:
        raise ValueError('Result is not a clean same-author review packet')
    inventory=(ROOT/'development/final_precard/evidence/TEST_INVENTORY.json').read_bytes()
    nodes=json.loads(inventory)
    if len(nodes)!=len(set(nodes)) or hashlib.sha256(inventory).hexdigest()!=m['test_inventory_sha256']:
        raise ValueError('Test inventory digest/uniqueness mismatch')
    return {'status':'READY FOR EXTERNAL CHECK','public_files_verified':len(m['public_files']),
            'work_branch_sha':m['work_branch_sha'],'results_sha256':m['results_sha256'],
            'private_kit_replay':'requires separately supplied pinned archive',
            'external_approval':False,'card_1_authorized':False}


if __name__=='__main__': print(json.dumps(verify(),sort_keys=True))
