#!/usr/bin/env python3
"""Replay the bounded packet in a new directory, preserving posted originals.

Usage: python reproduce.py /absolute/path/to/new-replay-directory
The destination must not already exist. No network or installation is used.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    target = destination/'research/reality_framework'/HERE.name
    target.parent.mkdir(parents=True)
    shutil.copytree(HERE, target, ignore=shutil.ignore_patterns('__pycache__','MANIFEST.json','REPLAY_COMPARISON.json','rg_results.json'))
    relative_profile = Path('grut_candidate_1_20261018/code/results_profile.npz')
    profile_target = target.parent/relative_profile
    profile_target.parent.mkdir(parents=True)
    shutil.copy2(HERE.parent/relative_profile, profile_target)
    commands = ['rg_symbolic.py','effective_potential.py','flow_check.py','verify.py','plot_results.py']
    execution = []
    for script in commands:
        start = time.perf_counter()
        proc = subprocess.run([sys.executable,str(target/script)],text=True,capture_output=True)
        (target/(script.removesuffix('.py')+'_replay.log')).write_text(proc.stdout+proc.stderr)
        execution.append({'command':[sys.executable,str(target/script)],'returncode':proc.returncode,'elapsed_seconds':time.perf_counter()-start})
        if proc.returncode:
            raise RuntimeError(f'{script} failed; see {target}')
    # Compare physical/algebraic outputs; timings, paths and platform strings are
    # execution metadata and are allowed to differ.
    checks = []
    def compare(filename, keys):
        old = json.loads((HERE/filename).read_text())
        new = json.loads((target/filename).read_text())
        for key in keys:
            a,b=old[key],new[key]
            if isinstance(a,(list,dict)):
                same = a == b
            else:
                same = bool(np.isclose(a,b,rtol=1e-12,atol=1e-15)) if isinstance(a,(float,int)) else a==b
            checks.append({'file':filename,'key':key,'passed':same})
    compare('rg_results.json',['checks_passed','full_identity_residual_monomials','beta_B_expressions','matching_beta'])
    compare('potential_results.json',['vacuum','bulk','classical_profile_flags','gravity_vacuum_shift'])
    compare('flow_results.json',['beta_at_matching','trajectories','scale_variations','bookkeeping','checks'])
    compare('VERIFICATION.json',['checks_passed','total_checks','all_checks_pass'])
    result={'source':str(HERE),'replay':str(target),'execution':execution,'checks':checks,
            'all_checks_pass':all(c['passed'] for c in checks),
            'profile_sha256':hashlib.sha256(profile_target.read_bytes()).hexdigest(),
            'limits':'Internal replay on the available numerical runtime, not an external replication or empirical validation.'}
    (target/'REPLAY_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'replay':str(target),'comparisons':len(checks),'all_checks_pass':result['all_checks_pass']}))
    if not result['all_checks_pass']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
