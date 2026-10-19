#!/usr/bin/env python3
"""Finite canonical checks for the rotating fixed-charge determinant formula."""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent


def clean(p):
    return {k: F(v) for k,v in p.items() if v}


def add(p,q):
    r=p.copy()
    for k,v in q.items(): r[k]=r.get(k,F(0))+v
    return clean(r)


def mul(p,q):
    r={}
    for a,x in p.items():
        for b,y in q.items(): r[a+b]=r.get(a+b,F(0))+x*y
    return clean(r)


def determinant_polynomial(K,G):
    n=len(K)
    entries=[[clean({0:K[i][j],1:G[i][j],2:F(int(i==j))}) for j in range(n)] for i in range(n)]
    result={}
    for perm in itertools.permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term={0:F(sign)}
        for i,j in enumerate(perm): term=mul(term,entries[i][j])
        result=add(result,term)
    return result


def generator(K,G):
    n=K.shape[0]
    return np.block([[np.zeros((n,n)),np.eye(n)],[-K,-G]])


def positive_nonzero(values):
    return np.sort(values.imag[(values.imag>1e-6)&(abs(values.real)<1e-7)])


def main():
    protocol_path=HERE/'theory_protocol.json'
    raw=protocol_path.read_bytes()
    protocol=json.loads(raw)
    output=HERE/'theory_toy_results.json'
    if output.exists(): raise SystemExit('Recorded toy output exists; refusing to overwrite.')
    O=F(4,5)
    k=F(-1,2)
    simple_nu2=k+4*O*O
    assert simple_nu2>0
    Kf=[[F(0) for _ in range(4)] for _ in range(4)]
    Kf[0][0]=k;Kf[0][1]=Kf[1][0]=F(1,10);Kf[1][1]=F(7,10);Kf[3][3]=F(2,5)
    Gf=[[F(0) for _ in range(4)] for _ in range(4)]
    for i in (0,1): Gf[i][i+2]=2*O;Gf[i+2][i]=-2*O
    keep=[0,1,3]
    Kr=[[Kf[i][j] for j in keep] for i in keep]
    Gr=[[Gf[i][j] for j in keep] for i in keep]
    Kr[0][0]+=4*O*O
    Kv=[[F(int(i==j))*(1-O*O) for j in range(4)] for i in range(4)]
    df=determinant_polynomial(Kf,Gf)
    dr=determinant_polynomial(Kr,Gr)
    dv=determinant_polynomial(Kv,Gf)
    identity=add(df,{p+2:-c for p,c in dr.items()})
    assert not identity,identity
    K=np.array(Kf,float);G=np.array(Gf,float)
    R=np.array(Kr,float);H=np.array(Gr,float)
    V=np.array(Kv,float)
    evR=np.linalg.eigvalsh(R)
    assert min(evR)>0
    Af=generator(K,G);Ar=generator(R,H);Av=generator(V,G)
    ef=np.linalg.eigvals(Af);er,vr=np.linalg.eig(Ar);ev=np.linalg.eigvals(Av)
    ff=positive_nonzero(ef);fr=positive_nonzero(er);fv=positive_nonzero(ev)
    assert len(ff)==len(fr)==3
    assert len(fv)==4
    assert max(abs(ff-fr))<1e-10
    assert max(abs(er.real))<1e-10
    reconstruct=[]
    for index,s in enumerate(er):
        if s.imag<=1e-6: continue
        qr=vr[:3,index]
        theta=2*float(O)*qr[0]/s
        q=np.array([qr[0],qr[1],theta,qr[2]])
        z=np.r_[q,s*q]
        residual=np.linalg.norm(Af@z-s*z)/np.linalg.norm(z)
        charge=2*float(O)*q[0]-s*q[2]
        assert residual<1e-10
        assert abs(charge)<1e-10
        reconstruct.append({'frequency':s.imag,'generator_residual':residual,'charge_abs':abs(charge)})
    rng=np.random.default_rng(20261019)
    energy_checks=[]
    for _ in range(4):
        qr=rng.normal(size=3);velocity=rng.normal(size=3);theta=rng.normal()
        q=np.array([qr[0],qr[1],theta,qr[2]])
        vel=np.array([velocity[0],velocity[1],2*float(O)*qr[0],velocity[2]])
        full=.5*(vel@vel+q@K@q)
        reduced=.5*(velocity@velocity+qr@R@qr)
        assert abs(full-reduced)<1e-12
        energy_checks.append({'full':full,'reduced':reduced,'absolute_difference':abs(full-reduced)})

    delta=add(df,{p:-c for p,c in dv.items()})
    def val(p,z):
        return sum(float(c)*z**power for power,c in p.items())
    def integrand(z):
        if z<1: return math.log(val(df,z))-math.log(val(dv,z))
        return math.log1p(val(delta,z)/val(dv,z))
    low,elow=quad(integrand,0,1,epsabs=2e-10,epsrel=2e-10,limit=200)
    high,ehigh=quad(integrand,1,np.inf,epsabs=2e-10,epsrel=2e-10,limit=200)
    integral=(low+high)/(2*math.pi)
    frequency_sum=.5*(sum(fr)-sum(fv))
    assert abs(integral-frequency_sum)<2e-9
    wrong_reduced_uv_degree=max(dr)-max(dv)
    result={
        'protocol_sha256':hashlib.sha256(raw).hexdigest(),
        'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
        'two_coordinate_toy':{'Omega_exact':str(O),'k_exact':str(k),'fixedQ_frequency_squared_exact':str(simple_nu2),
                              'unfixed_charge_curvature_exact':str(k/simple_nu2),
                              'meaning':'Negative grandcanonical charge curvature; fixed charge is a stable oscillator, so an unfixed real Gaussian trace is not the prescription.'},
        'full_determinant_polynomial':{str(p):str(c) for p,c in sorted(df.items())},
        'reduced_determinant_polynomial':{str(p):str(c) for p,c in sorted(dr.items())},
        'vacuum_determinant_polynomial':{str(p):str(c) for p,c in sorted(dv.items())},
        'exact_full_minus_zeta2_reduced_residual':{},
        'reduced_stiffness_eigenvalues':evR.tolist(),
        'full_nonzero_positive_frequencies':ff.tolist(),'reduced_positive_frequencies':fr.tolist(),
        'vacuum_positive_frequencies':fv.tolist(),
        'full_smallest_absolute_generator_roots':sorted(abs(ef).tolist())[:2],
        'max_full_reduced_frequency_difference':float(max(abs(ff-fr))),
        'reconstructed_mode_checks':reconstruct,'energy_checks':energy_checks,
        'logdet_integral_over_2pi':integral,'quadrature_error_estimate_over_2pi':(elow+ehigh)/(2*math.pi),
        'half_frequency_difference':frequency_sum,'integral_frequency_difference':integral-frequency_sum,
        'wrong_reduced_vs_full_vacuum_logdet_UV_power':wrong_reduced_uv_degree,
        'normalization_conclusion':'Multiply the charge-reduced determinant by zeta^2 before comparing to the original full vacuum determinant. No additional I or Qprime constant in an orthonormal kinetic basis. Phase zero energy and collective normalization are not a deletion of an arbitrary vacuum oscillator.',
        'limits':protocol['scope_limits'],
    }
    output.write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'theory_toy_execution.log').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
