#!/usr/bin/env python3
"""Read-only integrity, independent identities, and optional fresh replay."""
import argparse
import cmath
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np
from scipy.integrate import simpson

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--replay',action='store_true')
args=parser.parse_args()
R=json.loads((ROOT/'RESULTS.json').read_text())
I=json.loads((ROOT/'INPUTS.json').read_text())
manifest=json.loads((ROOT/'MANIFEST.json').read_text())['files']
actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
        and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
assert actual==set(manifest), ('manifest key set',actual^set(manifest))
for name,h in manifest.items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h, name
assert R['failed_checks']==0 and all(x['passed'] for x in R['checks'])
assert R['check_count']==len(R['checks'])==41

with (ROOT/'data/endpoint_map.csv').open() as f:
    rows=list(csv.DictReader(f))
assert len(rows)==len(R['histories'])==34
for a,b in zip(rows,R['histories']):
    for k,v in a.items():
        assert float(v)==b[k],(k,v,b[k])

raw_checks=[]
for file in sorted((ROOT/'data').glob('history_*.csv')):
    d=np.genfromtxt(file,delimiter=',',names=True)
    assert len(d)==I['raw_points']
    assert np.all(d['Omega_radiation']>0)
    assert np.max(abs(d['Omega_radiation']+d['Omega_kinetic']+d['Omega_potential']-1))<2e-14
    assert np.max(abs(d['du_dN']))<math.sqrt(6)/I['b']
    assert np.max(abs(d['u']-d['u'][0]))<=math.sqrt(6)/I['b']*d['N'][-1]+1e-12
    # Separate quadrature of Raychaudhuri, using only archived grid samples.
    log_ratio=math.log(d['H_GeV'][0]/d['H_GeV'][-1])
    integral=simpson(d['epsilon'],x=d['N'])
    energy_error=abs(integral-log_ratio)
    assert energy_error<2e-6,(file.name,energy_error)
    row={'file':file.name,'raw_energy_logH_error':float(energy_error)}
    # Integrating-factor Klein--Gordon identity, not the algebraic Friedmann constraint.
    if d['u'][0]>=60:
        u=d['u'];nn=d['N'];H=d['H_GeV']
        force=2*R['theory']['rho_star_GeV4']*np.exp(u)*np.expm1(u)
        P=np.exp(3*nn)*H*d['du_dN']
        Q=simpson(np.exp(3*nn)*force/(I['b']**2*I['M_GeV']**2*H),x=nn)
        error=abs((P[-1]-P[0]+Q)/max(abs(P[-1]),abs(Q)))
        assert error<2e-6,(file.name,error)
        row['raw_momentum_integral_relative_error']=float(error)
    raw_checks.append(row)

# Analytic eigenvalues independent of the calculation's finite-difference Jacobian.
for p in R['phase_plane']:
    lam2=p['lambda_exponential']**2
    if p['b']<1:
        root=cmath.sqrt(64/lam2-15)
        analytic=[-.5*(1+root),-.5*(1-root)]
        assert abs(p['scalar_fraction']-p['b']**2)<1e-14
    else:
        analytic=[complex((lam2-6)/2),complex(lam2-4)]
        assert abs(p['total_w']-(lam2/3-1))<1e-14
    numerical=[complex(*x) for x in p['jacobian_eigenvalues']]
    assert min(max(abs(numerical[0]-analytic[0]),abs(numerical[1]-analytic[1])),
               max(abs(numerical[0]-analytic[1]),abs(numerical[1]-analytic[0])))<2e-8

comparison={'manifest_files':len(manifest),'internal_checks':R['check_count'],'raw_identities':raw_checks,
            'analytic_phase_eigenvalues':'pass'}


def scientific_leaves(obj,path=''):
    if isinstance(obj,dict):
        out={}
        for k,v in obj.items():
            if k=='evaluation_count':
                continue
            out.update(scientific_leaves(v,path+'/'+k))
        return out
    if isinstance(obj,list):
        out={}
        for i,v in enumerate(obj):
            out.update(scientific_leaves(v,path+'/'+str(i)))
        return out
    return {path:obj}


if args.replay:
    with tempfile.TemporaryDirectory(prefix='grut-selection-replay-') as td:
        dest=Path(td)
        shutil.copytree(ROOT/'code',dest/'code',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(ROOT/'INPUTS.json',dest/'INPUTS.json')
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
        proc=subprocess.run([sys.executable,str(dest/'code/calculate.py')],capture_output=True,text=True,env=env)
        assert proc.returncode==0,proc.stdout+proc.stderr
        a=scientific_leaves(R)
        fresh=json.loads((dest/'RESULTS.json').read_text())
        z=scientific_leaves(fresh)
        assert a.keys()==z.keys()
        max_relative=0.
        for key,value in a.items():
            new=z[key]
            if isinstance(value,(int,float)) and not isinstance(value,bool):
                assert math.isclose(value,new,rel_tol=3e-8,abs_tol=1e-100),(key,value,new)
                max_relative=max(max_relative,abs(value-new)/max(abs(value),abs(new),1e-300))
            else:
                assert value==new,(key,value,new)
        oldfiles={p.name for p in (ROOT/'data').glob('*.csv')}
        newfiles={p.name for p in (dest/'data').glob('*.csv')}
        assert oldfiles==newfiles
        for file in oldfiles:
            old=np.genfromtxt(ROOT/'data'/file,delimiter=',',names=True)
            new=np.genfromtxt(dest/'data'/file,delimiter=',',names=True)
            assert old.dtype.names==new.dtype.names and old.shape==new.shape
            for col in old.dtype.names:
                assert np.allclose(old[col],new[col],rtol=3e-8,atol=1e-100),(file,col)
        comparison.update(replay_scientific_leaves=len(a),replay_max_relative_difference=max_relative,
                          replay_raw_arrays=len(oldfiles))
print(json.dumps(comparison,indent=2))
