"""Packet integrity, analytic charge floor and declared numerical gates.

These are internal checks; they are not independent scientific replication.
"""
from pathlib import Path
import hashlib,json,tarfile
import numpy as np
HERE=Path(__file__).resolve().parent

def main():
    checks=[]
    def check(name,ok,**v):checks.append({'name':name,'pass':bool(ok),**v})
    f=json.loads((HERE/'SOURCE_FREEZE.json').read_text())
    for p,meta in f['files'].items():check('frozen_'+p,hashlib.sha256((HERE/p).read_bytes()).hexdigest()==meta['sha256'])
    gf=json.loads((HERE/'GAUSSIAN_SOURCE_FREEZE.json').read_text())
    for p,sha in gf['files'].items():check('gaussian_frozen_'+p,hashlib.sha256((HERE/p).read_bytes()).hexdigest()==sha)
    v=json.loads((HERE/'VALIDATION.json').read_text());gv=json.loads((HERE/'GAUSSIAN_VALIDATION.json').read_text())
    check('current_solver_validation',v['failed']==0,passed=v['passed'])
    check('gaussian_wrapper_validation',gv['failed']==0,passed=gv['passed'])
    for name,expected in [('FAILED_VALIDATION_01.json',2),('FAILED_VALIDATION_02.json',1)]:
        x=json.loads((HERE/name).read_text());check('retained_'+name,x['failed']==expected,failed=x['failed'])
    r=json.loads((HERE/'RESULTS.json').read_text())
    check('all_declared_main_jobs_disposed',len([x for x in r['status'] if not x['name'].startswith('gaussian')])==11)
    check('all_declared_gaussian_jobs_disposed',len([x for x in r['status'] if x['name'].startswith('gaussian')])==3)
    for name,x in r['runs'].items():
        check('numerical_gate_'+name,x['numerical_gates_pass'] and x['complete'],charge_drift=x['charge_relative_drift'],work_residual=x['conformal_budget_residual'])
        d=np.genfromtxt(HERE/'raw'/(name+'.csv'),delimiter=',',names=True)
        # rho has unit lambda_phi v^4; nbar=n_Q/(sqrt(lambda_phi) v^3).
        # m*/m_X=(lambda_X/lambda_phi)^(1/4)=0.1^(1/4).
        floor=(.1**.25)*np.abs(d['Q_com'])/(x['L']**3*d['a']**3)
        margin=float(np.min(d['rho']-floor))
        check('classical_charge_floor_'+name,margin>=-1e-12,minimum_dimensionless_margin=margin)
    raw=json.loads((HERE/'RAW_MANIFEST.json').read_text())
    with tarfile.open(HERE/raw['archive'],'r:gz') as tf:
        members={m.name:m for m in tf.getmembers() if m.isfile()}
        expected={x['path']:x for x in raw['files'] if x['archived']}
        check('strict_archive_leaf_paths',set(members)==set(expected))
        for path,entry in expected.items():
            payload=tf.extractfile(members[path]).read()
            check('archive_'+path,len(payload)==entry['bytes'] and hashlib.sha256(payload).hexdigest()==entry['sha256'])
    if (HERE/'MANIFEST.json').exists():
        entries=json.loads((HERE/'MANIFEST.json').read_text())
        for e in entries['files']:check('manifest_'+e['path'],hashlib.sha256((HERE/e['path']).read_bytes()).hexdigest()==e['sha256'])
    out={'checks':checks,'passed':sum(c['pass'] for c in checks),'failed':sum(not c['pass'] for c in checks),
         'scope':'Integrity and internal mathematical/numerical checks; no outside replication and no empirical validation.'}
    (HERE/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':out['passed'],'failed':out['failed'],'failures':[c for c in checks if not c['pass']]},indent=2))
    if out['failed']:raise SystemExit(1)

if __name__=='__main__':main()
