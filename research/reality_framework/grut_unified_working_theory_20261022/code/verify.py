#!/usr/bin/env python3
"""Read-only raw identities, Numerov control, quantum readout and fresh replay."""
import argparse
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
from scipy.interpolate import CubicSpline
from scipy.stats import binom

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');args=parser.parse_args()
I=json.loads((ROOT/'INPUTS.json').read_text());R=json.loads((ROOT/'RESULTS.json').read_text())
files=json.loads((ROOT/'MANIFEST.json').read_text())['files']
actual={str(x.relative_to(ROOT)) for x in ROOT.rglob('*') if x.is_file() and x.name!='MANIFEST.json' and '__pycache__' not in x.parts}
assert set(files)==actual
for f,sha in files.items():
    assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==sha,f
source=ROOT.parent/I['source_profile_relative_to_parent']
assert hashlib.sha256(source.read_bytes()).hexdigest()==I['source_profile_sha256']
assert [I[k] for k in ['g_X','lambda_phi','lambda_X','h_D','m_D0_over_gv','lambda_D','v_over_reduced_M']]==[.1,.01,.001,.01,1.,.001,1e-5]
assert R['failed_checks']==0 and all(x['passed'] for x in R['checks'])
assert len(R['checks'])==R['check_count']==51
stats={x['Q']:x for x in R['profiles']}
raw_checks=[];phase_data={}


def numerov(spline,L,k,h):
    """Independent finite-difference scattering from the archived profile."""
    n=round((L+10)/h);rr=np.arange(n+1)*h
    ss=np.where(rr<=L,spline(np.minimum(rr,L)),1.)
    q=k*k-(I['h_D']/I['g_X'])**2*(ss*ss-1)
    u=np.zeros(n+1);u[1]=h
    for j in range(1,n):
        u[j+1]=(2*(1-5*h*h*q[j]/12)*u[j]-(1+h*h*q[j-1]/12)*u[j-1])/(1+h*h*q[j+1]/12)
    j=n-2
    up=(u[j-2]-8*u[j-1]+8*u[j+1]-u[j+2])/(12*h)
    return (math.atan2(k*u[j],up)-k*rr[j]+math.pi/2)%math.pi-math.pi/2


for Q in I['main_charges']:
    d=np.genfromtxt(ROOT/f'data/profile_Q{Q}.csv',delimiter=',',names=True)
    r=d['rho'];F=d['F'];S=d['S'];om=stats[Q]['Omega']
    aa=I['lambda_phi']/I['g_X']**2;bb=I['lambda_X']/I['g_X']**2
    potential=aa/4*(S*S-1)**2+.5*S*S*F*F+bb/4*F**4
    density=.5*(d['Fprime']**2+d['Sprime']**2)+.5*om*om*F*F+potential
    E=4*math.pi*simpson(r*r*density,x=r)
    charge=om*4*math.pi*simpson(r*r*F*F,x=r)/I['g_X']**2
    depletion=aa*4*math.pi*simpson(r*r*(1-S*S),x=r)
    assert abs(E/stats[Q]['Ebar']-1)<2e-8
    assert abs(charge/Q-1)<2e-8
    assert abs(E/depletion-1)<2e-8
    assert stats[Q]['exact_scattering_length_lower_bound']<=stats[Q]['scattering_length_over_inverse_gv']<=stats[Q]['exact_scattering_length_upper_bound']
    curves=np.genfromtxt(ROOT/f'data/phases_Q{Q}.csv',delimiter=',',names=True)
    assert len(curves)==I['probe_k_points']
    phase_data[Q]=curves
    spline=CubicSpline(r,S)
    errors=[]
    for index in [0,len(curves)//2,len(curves)-1]:
        k=float(curves['k_over_gv'][index]);target=float(curves['delta_l0'][index])
        coarse=numerov(spline,r[-1],k,.02);fine=numerov(spline,r[-1],k,.01)
        err=abs(fine-target)
        assert err<2e-7,(Q,k,coarse,fine,target)
        assert abs(coarse-fine)<2e-7
        errors.append(err)
    raw_checks.append(dict(Q=Q,raw_mass_relative_error=float(abs(E/stats[Q]['Ebar']-1)),
                           raw_sumrule_relative_error=float(abs(E/depletion-1)),
                           max_Numerov_s_wave_phase_error=max(errors)))

# Finite mode assembly uses explicitly normalized Simpson weights; no analytic
# two-state distance formula is used for this first single-probe check.
Q0,Q1=I['readout_charges'];d0=phase_data[Q0];d1=phase_data[Q1];ks=d0['k_over_gv']
coeff=np.ones(len(ks));coeff[1:-1:2]=4;coeff[2:-1:2]=2
w=coeff*(ks[1]-ks[0])/3*np.exp(-.5*((ks-I['packet_k_center'])/I['packet_k_sigma'])**2)
w/=w.sum()
psi=[]
for d in [d0,d1]:
    a=np.stack([np.exp(2j*d['delta_l0']),np.exp(2j*d['delta_l1'])],axis=1)*np.sqrt(w[:,None]/2)
    psi.append(a.ravel())
C=np.vdot(psi[0],psi[1]);rho0=np.outer(psi[0],psi[0].conj());rho1=np.outer(psi[1],psi[1].conj())
distance=float(np.abs(np.linalg.eigvalsh(rho0-rho1)).sum()/2)
assert abs((1+distance)/2-R['readout']['single_probe_optimal_success'])<2e-12
assert abs(abs(C)**2-R['readout']['overlap_squared'])<2e-14
# One monochromatic partial wave supplies only a global phase.
one0=np.array([np.exp(2j*d0['delta_l0'][30])]);one1=np.array([np.exp(2j*d1['delta_l0'][30])])
assert np.max(abs(np.outer(one0,one0.conj())-np.outer(one1,one1.conj())))<5e-16

p0,p1=R['readout']['measurement_plus_probabilities']
N=R['readout']['fixed_measurement_probe_count_95']
for n in [1,10,N-1,N]:
    kk=np.arange(n+1)
    success=float(.5*np.maximum(np.exp(binom.logpmf(kk,n,p0)),np.exp(binom.logpmf(kk,n,p1))).sum())
    if n==N:
        assert abs(success-R['readout']['fixed_measurement_success_at_N'])<2e-10 and success>=.95
    if n==N-1:
        assert abs(success-R['readout']['fixed_measurement_success_before_N'])<2e-10 and success<.95

out=dict(manifest_files=len(files),internal_checks=R['check_count'],raw_profile_and_Numerov_controls=raw_checks,
         independently_assembled_single_probe_success=(1+distance)/2,explicit_binomial_enumeration='pass')


def leaves(obj,p=''):
    if isinstance(obj,dict):
        result={}
        for k,v in obj.items():
            if k!='evaluation_count':result.update(leaves(v,p+'/'+k))
        return result
    if isinstance(obj,list):
        result={}
        for j,v in enumerate(obj):result.update(leaves(v,p+'/'+str(j)))
        return result
    return {p:obj}


if args.replay:
    with tempfile.TemporaryDirectory(prefix='grut-common-model-') as td:
        dest=Path(td)/'packet';dest.mkdir()
        shutil.copytree(ROOT/'code',dest/'code',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(ROOT/'INPUTS.json',dest/'INPUTS.json')
        sp=Path(td)/I['source_profile_relative_to_parent'];sp.parent.mkdir(parents=True)
        shutil.copy2(source,sp)
        proc=subprocess.run([sys.executable,str(dest/'code/calculate.py')],capture_output=True,text=True,
                             env=dict(os.environ,OPENBLAS_NUM_THREADS='1'))
        assert proc.returncode==0,proc.stdout+proc.stderr
        fresh=json.loads((dest/'RESULTS.json').read_text());a=leaves(R);b=leaves(fresh)
        assert a.keys()==b.keys()
        worst=0.
        for p,x in a.items():
            y=b[p]
            if isinstance(x,(float,int)) and not isinstance(x,bool):
                assert math.isclose(x,y,rel_tol=3e-7,abs_tol=2e-12),(p,x,y)
                worst=max(worst,abs(x-y)/max(abs(x),abs(y),1e-300))
            else:assert x==y,(p,x,y)
        names={x.name for x in (ROOT/'data').glob('*.csv')}
        assert names=={x.name for x in (dest/'data').glob('*.csv')}
        for name in names:
            a0=np.genfromtxt(ROOT/'data'/name,delimiter=',',names=True)
            a1=np.genfromtxt(dest/'data'/name,delimiter=',',names=True)
            assert a0.shape==a1.shape and a0.dtype.names==a1.dtype.names
            for col in a0.dtype.names:assert np.allclose(a0[col],a1[col],rtol=3e-7,atol=2e-12),(name,col)
        out.update(replay_scientific_leaves=len(a),replay_worst_relative_difference=worst,replay_raw_arrays=len(names))
print(json.dumps(out,indent=2))
