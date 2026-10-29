"""Exact piecewise-linear radial integration of the parent bin-center model.

Each speed has its own thin shell. Sort every retarded-time knot and sweep
piecewise-linear source amplitudes so radial integration resolves thin shells
without taking a fixed radial step. This is a sensitivity model, not a claim
that the physical source consists of spectral delta functions.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
from scipy.integrate import simpson

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'grut_m1_local_cloud_20261029'
sys.path.insert(0,str(OLD/'code'))
from local_rates import pulse_series
NC=1.2141071422639381e-15

def sweep(T,t,e,flux,r0):
    v=np.sqrt(e*e-1)/e;dt=t[1]-t[0];nt=len(t)
    x=r0+v[:,None]*(T-t[None,:]);F=flux.T
    slopes=-np.diff(F,axis=1)/(dt*v[:,None]**2)
    change=np.zeros_like(F);jump=np.zeros_like(F)
    change[:,-1]=slopes[:,-1];change[:,1:-1]=slopes[:,:-1]-slopes[:,1:];change[:,0]=-slopes[:,0]
    jump[:,-1]=F[:,-1]/v;jump[:,0]=-F[:,0]/v
    order=np.argsort(x.ravel());r=x.ravel()[order];change=change.ravel()[order];jump=jump.ravel()[order]
    m=np.cumsum(change);dx=np.diff(r)
    inc=np.r_[0.,m[:-1]*dx]+jump
    s=np.cumsum(inc)
    right=s[:-1]+m[:-1]*dx
    N=float(np.sum((s[:-1]+right)*.5*dx))
    # 3-point Gauss quadrature on each linear segment, with explicit r^-2.
    nodes=np.array([-.7745966692414834,0.,.7745966692414834]);weights=np.array([5/9,8/9,5/9])
    V=0.
    for q,w in zip(nodes,weights):
        f=(q+1)/2;sm=s[:-1]+f*m[:-1]*dx;rr=r[:-1]+f*dx
        V+=float(np.sum(w*.5*dx*sm*sm/(rr*rr)))/(4*np.pi)
    target=float(simpson(flux.sum(axis=1),x=t))
    cv=np.sqrt(NC*V)/(NC*N)
    return {'N_radial':N,'N_source_simpson':target,'relative_number_error':N/target-1,'density_cv':float(cv),
            'after_last_knot_residual':float(s[-1]),'min_segment_start':float(np.min(s)),'max_segment_start':float(np.max(s))}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='source_t400');p.add_argument('--dt',type=float,default=10.)
    p.add_argument('--tag',required=True);a=p.parse_args();z=np.load(OLD/(a.source+'.npz'))
    t,e,plus,minus=pulse_series(z,a.dt);de=float(np.diff(z['energy_edges'])[0]);plus*=de;r0=float(z['probe_radius_mX_inverse'])
    rows=[]
    for T in (100000.,200000.,500000.,1000000.,10000000.,100000000.,1000000000.,10000000000.):
        rows.append({'mX_t':T,'plus':sweep(T,t,e,plus,r0),'minus':sweep(T,t,e,minus,r0)})
    out={'configuration':vars(a),'rows':rows,
         'scope':'Exact knot-swept density moments for the same bin-center, delta-speed model used in the parent Poisson study. No claim that its discrete spectral lines are physical.'}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
