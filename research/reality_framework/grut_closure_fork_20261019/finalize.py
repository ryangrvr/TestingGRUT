#!/usr/bin/env python3
"""Packet integrity and recorded-result checks, not a novelty certificate."""
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    result=json.loads((HERE/'results.json').read_text())
    initial=json.loads((HERE/'results_initial.json').read_text())
    assert result['all_checks_pass']
    assert not initial['checks']['midpoint_residual_scales_quadratically']
    assert len(result['classical_BVP_runs'])==6
    for name,key in [('protocol.json','protocol_sha256'),('refinement_protocol.json','refinement_protocol_sha256'),('curvature_transfer_protocol.json','curvature_transfer_protocol_sha256')]:assert digest(HERE/name)==result[key]
    inputs=[HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz',HERE.parent/'grut_candidate_1_20261018/code/results.json',HERE.parent.parent.parent/'GRUT_III/calc/v11_quasi_ds_kms_defect.py',HERE.parent.parent.parent/'GRUT_III/calc/v14_massless_sector_first_order.py']
    assert digest(inputs[0])==result['profile_sha256']=='0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3'
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json')
    for p in files:
        if p.suffix=='.py':compile(p.read_text(),str(p),'exec')
        if p.suffix=='.json':json.loads(p.read_text())
        if p.suffix=='.md':
            for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
                if not target.startswith(('https://','http://','#')):assert (p.parent/target.split('#')[0]).exists(),(p.name,target)
    manifest={'parent_commit':'8fdd674d5e487dc050d685890e7c6defc35e18b1','branch':'codex/grut-closure-fork-20261019','files':[{'path':p.name,'bytes':p.stat().st_size,'sha256':digest(p)} for p in files],'source_inputs':[{'path':str(p.relative_to(HERE.parent.parent.parent)),'sha256':digest(p)} for p in inputs],'executed_conditions':result['checks'],'scope':'Byte integrity, syntax, links and recorded numerical outcomes. No independent external validation, cosmological path integral, new empirical fit or priority certificate.'}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'files':len(files),'classical_solves':6,'numerical_conditions_pass':all(result['checks'].values()),'initial_coarse_failure_preserved':True}))


if __name__=='__main__':main()
