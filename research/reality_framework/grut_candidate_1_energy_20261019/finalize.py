#!/usr/bin/env python3
"""Integrity/format checks and generated packet manifest; not a physics gate."""
import csv
import hashlib
import json
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json')
    checks=[]
    for p in files:
        if p.suffix=='.py':compile(p.read_text(),str(p),'exec');checks.append('Python syntax: '+p.name)
        if p.suffix=='.json':json.loads(p.read_text());checks.append('JSON parse: '+p.name)
        if p.suffix=='.jsonl':
            for line in p.read_text().splitlines():json.loads(line)
            checks.append('JSONL parse: '+p.name)
        if p.suffix=='.md':
            for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
                if not target.startswith(('https://','http://','#')):
                    assert (p.parent/target.split('#')[0]).exists(),(p.name,target)
            checks.append('Markdown relative links: '+p.name)
    with (HERE/'spectral_points.csv').open() as f:rows=list(csv.DictReader(f))
    result=json.loads((HERE/'spectral_results.json').read_text())
    assert len(rows)==len(result['points'])==2362
    verification=json.loads((HERE/'VERIFICATION.json').read_text()); assert verification['all_checks_pass']
    binding=json.loads((HERE/'binding_results.json').read_text()); assert binding['all_checks_pass']
    # The immutable parent profile/pole inputs are outside this packet but
    # explicitly part of its reproducibility dependency closure.
    inputs=[HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz',HERE.parent/'grut_candidate_1_20261018/code/results.json',HERE.parent/'grut_candidate_1_fluctuations_20261019/pole_results.json']
    manifest={'parent_commit':'1f370987f2fa96a52426ec95680a4f3ad7cbd373','intended_branch':'codex/grut-candidate-1-quantum-energy-20261019','files':[{'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in files],
      'immutable_parent_inputs':[{'path':str(p.relative_to(HERE.parent)),'sha256':sha(p)} for p in inputs],
      'format_checks':checks,'spectral_csv_rows':len(rows),'physics_verification_checks':verification['total_checks'],
      'executions':{'original_spectral_seconds':result['elapsed_seconds'],'strict_full_replay_seconds':json.loads((HERE/'audit_results.json').read_text())['elapsed_seconds'],'graphs_seconds':json.loads((HERE/'graphs_results.json').read_text())['elapsed_seconds'],'integer_integration_seconds':binding['execution_seconds']},
      'scope':'Byte integrity, syntax and link checks; not physical validation. No self-hash, ignored cache files or unrelated preexisting files included.'}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'files':len(files),'format_checks':len(checks),'spectral_rows':len(rows),'physics_checks':verification['total_checks'],'manifest_written':True}))


if __name__=='__main__':main()
