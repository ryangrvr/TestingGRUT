#!/usr/bin/env python3
"""Replay in a new directory without overwriting the posted packet.

Usage: python reproduce.py /absolute/path/to/new-replay-directory
"""
from pathlib import Path
import argparse,os,json,shutil,subprocess,sys,time
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('destination',type=Path)
    dest=parser.parse_args().destination.resolve();dest.mkdir(parents=True,exist_ok=False)
    target=dest/'research/reality_framework'/HERE.name;target.parent.mkdir(parents=True)
    shutil.copytree(HERE,target,ignore=shutil.ignore_patterns('__pycache__','MANIFEST.json','REPRODUCTION.json','REPLAY_COMPARISON.json','mode_initial_results.json','mode_initial_execution.jsonl'))
    parents=['grut_candidate_1_20261018/code/results_profile.npz','grut_candidate_1_quantum_20261018/flow_check.py','grut_candidate_1_quantum_20261018/effective_potential.py','grut_candidate_1_quantum_20261018/potential_results.json']
    for path in parents:
        out=target.parent/path;out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(HERE.parent/path,out)
    commands=[['modes.py'],['modes.py','--box-check'],['angular_bound.py'],['modes.py','--ell5-check'],['pole.py'],['neutral_poles.py'],['verify.py'],['plot_results.py']]
    env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    execution=[]
    for i,cmd in enumerate(commands):
        before=time.perf_counter();run=subprocess.run([sys.executable,str(target/cmd[0]),*cmd[1:]],env=env,capture_output=True,text=True,timeout=600)
        (target/f'replay_{i+1:02d}.log').write_text(run.stdout+run.stderr)
        execution.append({'command':[cmd[0],*cmd[1:]],'returncode':run.returncode,'elapsed_seconds':time.perf_counter()-before})
        if run.returncode:raise RuntimeError(f'{cmd} failed; see {target}')
    checks=[]
    def check(name,value):checks.append({'name':name,'passed':bool(value)})
    old=json.loads((HERE/'mode_results.json').read_text());new=json.loads((target/'mode_results.json').read_text())
    check('40 spectra reproduced',len(old['spectra'])==len(new['spectra'])==40)
    check('no budget or eigensolver failure in replay',not new['skipped'] and all(t.get('status','converged')=='converged' for t in new['attempts']))
    key=lambda s:(s['method'],s['N'],s['L'],s['ell'])
    bykey={key(s):s for s in new['spectra']}
    for s in old['spectra']:
        t=bykey[key(s)]
        check('energy spectra '+str(key(s)),np.allclose(s['energy_checks']['real_amplitude_lowest_energy_eigenvalues'],t['energy_checks']['real_amplitude_lowest_energy_eigenvalues'],rtol=1e-7,atol=1e-9))
        oldfreq=[r['s_imag'] for r in s['smallest_roots'] if r['symmetry'] is None and r['s_imag']>1e-5]
        newfreq=[r['s_imag'] for r in t['smallest_roots'] if r['symmetry'] is None and r['s_imag']>1e-5]
        check('lowest nonsymmetry frequency '+str(key(s)),bool(oldfreq and newfreq) and np.isclose(min(oldfreq),min(newfreq),rtol=1e-7,atol=1e-9))
    for filename,keys in [('pole_results.json',['strict_pole_squared_over_v2','strict_pole_mass_over_v','strict_curvature_squared_over_v2']),('neutral_pole_results.json',[])]:
        a=json.loads((HERE/filename).read_text());b=json.loads((target/filename).read_text())
        if keys:
            for k in keys:check(filename+' '+k,np.isclose(a['matching_result'][k],b['matching_result'][k],rtol=1e-12,atol=1e-14))
        else:
            check('all four tensor pole masses',np.allclose([r['strict_pole_squared_over_v2'] for r in a['results']],[r['strict_pole_squared_over_v2'] for r in b['results']],rtol=1e-12,atol=1e-14))
        check(filename+' acceptance checks',b['all_checks_pass'])
    a=json.loads((HERE/'angular_bound_results.json').read_text());b=json.loads((target/'angular_bound_results.json').read_text())
    check('angular bound diagnostic reproduced',np.isclose(a['optimized']['optimized_B'],b['optimized']['optimized_B'],rtol=1e-10,atol=1e-10) and b['half_line_bound_certified'] is False)
    v=json.loads((target/'VERIFICATION.json').read_text());check('all replay cross-checks pass',v['all_checks_pass'])
    result={'source':str(HERE),'replay':str(target),'execution':execution,'checks':checks,'all_checks_pass':all(c['passed'] for c in checks),
            'scope':'Same-runtime numerical replay with tolerances, not bitwise identity of sparse eigenvectors or external validation.'}
    (target/'REPLAY_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'replay':str(target),'comparisons':len(checks),'all_checks_pass':result['all_checks_pass']}))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
