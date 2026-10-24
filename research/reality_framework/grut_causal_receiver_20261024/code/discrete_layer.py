#!/usr/bin/env python3
"""Exact incoming/evanescent lead matching for the implemented staggered grid."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

ROOT=Path(__file__).resolve().parents[1]
A0=-365/144;A1=87/64;A2=-3/32;A3=1/576


def lattice(k,h):
    sine=np.sin(k*h/2)
    return 2*sine/h*(1+sine*sine/6)


def reflection(k,h=.25,width=90.,sigma=20.,alpha=.1,mass=1.):
    nlayer=round(width/h);lead=5;n=nlayer+lead;width=nlayer*h
    momentum=lattice(k,h);omega=np.hypot(mass,momentum);s=-1j*omega
    roots=np.roots([A3,A2,A1,A0+(momentum*h)**2,A1,A2,A3])
    evanescent=sorted([z for z in roots if abs(z)>1.01],key=abs)
    assert len(evanescent)==2,(k,h,roots)
    zi=np.exp(1j*k*h);zr=zi.conjugate()
    matrix=lil_matrix((n+3,n+3),dtype=complex);rhs=np.zeros(n+3,dtype=complex)
    def coefficient(x):
        strength=sigma*np.clip(x/width,0,1)**3
        return (s+alpha)/(s+alpha+strength)
    def sample(j):
        if j==n:return {},0j
        if j>n:return {2*n-j:-1.},0j
        if j>=0:return {j:1.},0j
        return {n:zr**j,n+1:evanescent[0]**j,n+2:evanescent[1]**j},zi**j
    for j in range(n):
        cn=coefficient((j-lead)*h)
        for face,d in [(j,9/8),(j-1,-9/8),(j+1,-1/24),(j-2,1/24)]:
            folded=face if face<n else 2*n-1-face
            cf=coefficient((folded+.5-lead)*h)
            for node,g in [(face+1,9/8),(face,-9/8),(face+2,-1/24),(face-1,1/24)]:
                value=cn*d*cf*g/h**2;weights,constant=sample(node)
                for col,w in weights.items():matrix[j,col]+=value*w
                rhs[j]-=value*constant
        matrix[j,j]+=momentum**2
    for j in range(3):
        matrix[n+j,j]=1;matrix[n+j,n]=-zr**j
        matrix[n+j,n+1]=-evanescent[0]**j;matrix[n+j,n+2]=-evanescent[1]**j
        rhs[n+j]=zi**j
    answer=spsolve(matrix.tocsr(),rhs);r=answer[n]
    ideal=np.exp(-sigma*width*momentum*omega/(alpha**2+omega**2))
    return dict(k=float(k),lattice_momentum=float(momentum),mass=float(mass),h=h,width=width,
        sigma_max=sigma,alpha=alpha,reflected_energy=float(abs(r)**2),
        reflection_real=float(r.real),reflection_imag=float(r.imag),
        ideal_continuum_wall_energy=float(ideal),linear_system_residual=float(np.linalg.norm(matrix@answer-rhs,np.inf)))


def main():
    rows=[];controls=[]
    for h in [.5,.25,.125,.0625]:
        for mass in [1.,np.sqrt(6.)]:
            for k in [.01,.02,.03,.05,.1,.2,.5,1.,2.,3.]:
                rows.append(reflection(k,h,mass=mass))
        for k in [.01,.1,1.,3.]:
            r=reflection(k,h,sigma=0.)
            n=round(r['width']/h)+5;exact=-np.exp(2j*k*h*n)
            r['zero_layer_complex_error']=float(abs(complex(r['reflection_real'],r['reflection_imag'])-exact))
            controls.append(r)
    out=dict(scope='Actual staggered fourth-order spatial stencil, causal CFS frequency response, outer wall, exact oscillatory and evanescent lattice leads. No continuum-only entrance approximation.',
        rows=rows,zero_layer_controls=controls,
        maximum_sampled_reflection=max(r['reflected_energy'] for r in rows),
        maximum_sampled_reflection_charged=max(r['reflected_energy'] for r in rows if r['mass']==1.),
        one_percent_sampled_gate_passed=all(r['reflected_energy']<.01 for r in rows))
    (ROOT/'DISCRETE_LAYER_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('maximum sampled R',out['maximum_sampled_reflection'],'gate',out['one_percent_sampled_gate_passed'],flush=True)
    print('zero-layer maximum error',max(r['zero_layer_complex_error'] for r in controls),flush=True)
    for h in [.5,.25,.125,.0625]:
        for mass in [1.,np.sqrt(6.)]:
            group=[r for r in rows if r['h']==h and r['mass']==mass]
            print('h',h,'mass',mass,'max R',max(r['reflected_energy'] for r in group),flush=True)


def cold_controls():
    rows=[reflection(k,2.) for k in [.005,.01,.02,.03,.05,.1,.2]]
    out=dict(scope='The coarser spatial grid used only for the long soft free-packet benchmark; it is not a formation grid.',rows=rows)
    (ROOT/'COLD_GRID_CONTROLS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    if '--cold-grid' not in sys.argv:main()
    cold_controls()
