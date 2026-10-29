"""Collisionless Poisson cloud with a continuous, piecewise-linear E source.

Integrate the retarded source over emission time at each radius instead of
representing each energy bin by one speed. This resolves the spurious comb
of delta-speed shells when t becomes much larger than the source duration.
The flat threshold endcap is an assumed model, not measured physics.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
from scipy.integrate import simpson,quad

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'grut_m1_local_cloud_20261029'
sys.path.insert(0,str(OLD/'code'))
from local_rates import pulse_series
NC=1.2141071422639381e-15

def source_at_flux_energy(F,E,e,edges):
    de=e[1]-e[0];x=(E-e[0])/de
    j=np.floor(x).astype(np.int64);w=x-j
    w=np.where(x<=0,0,np.where(x>=len(e)-1,1,w));j=np.clip(j,0,len(e)-2)
    ti=np.arange(F.shape[0])[:,None]
    values=(1-w)*F[ti,j]+w*F[ti,j+1]
    return np.where((E>=edges[0])&(E<=edges[-1]),values,0.)

def profiles(T,t,e,edges,fluxes,r0,nv,chunk):
    if T<=t[-1]:raise ValueError('Observation time must exceed source duration')
    vmax=np.sqrt(1-1/edges[-1]**2)
    v0=np.linspace(0,vmax,nv);radius=r0+T*v0;S=np.empty((len(fluxes),nv))
    u=t[:,None];factor=1-u/T
    for i in range(0,nv,chunk):
        j=min(i+chunk,nv);v=v0[None,i:j]/factor
        valid=(v>=0)&(v<vmax)
        vv=np.clip(v,0,vmax);gamma=1/np.sqrt(1-vv*vv);E=gamma
        # E changes with u at fixed r: (dE/du)/v = gamma^3 v/(T-u).
        jac=np.where(valid,gamma**3*vv/(T-u),0.)
        for k,F in enumerate(fluxes):
            values=source_at_flux_energy(F,E,e,edges)
            S[k,i:j]=simpson(values*jac,x=t,axis=0)
    return radius,S

def calc(z,T,nv,dt,chunk):
    t,e,plus,minus=pulse_series(z,dt);edges=z['energy_edges'];de=edges[1]-edges[0]
    minus=minus/de;r0=float(z['probe_radius_mX_inverse'])
    r,S=profiles(T,t,e,edges,[plus,minus],r0,nv,chunk)
    N=np.array([simpson(a,x=r) for a in S]);
    # Piecewise-linear energy interpolation integrates by the trapezoid rule,
    # not Simpson's rule; the endpoint plateaus supply the two half bins.
    expect=np.array([np.trapezoid(simpson(a,x=t,axis=0),x=e) for a in (plus,minus)])
    expect+=.5*de*np.array([simpson(a[:,0]+a[:,-1],x=t) for a in (plus,minus)])
    pair=np.array([[simpson(S[i]*S[j]/(r*r),x=r)/(4*np.pi) for j in range(2)] for i in range(2)])
    cv=np.sqrt(NC*np.diag(pair))/(NC*N)
    ratio=1+NC*pair[0,1]/(NC*N[0]*NC*N[1])
    return {'mX_t':T,'N_via_radial_integral':N.tolist(),'N_via_source':expect.tolist(),
            'relative_number_error':(N/expect-1).tolist(),'density_cv':cv.tolist(),
            'mean_density_product_ratio':float(ratio),'max_velocity':float(np.sqrt(1-1/edges[-1]**2))}

def mono_shell(T,vel,duration,r0):
    a=r0+vel*(T-duration);b=r0+vel*T;width=vel*duration
    effective=4*np.pi*NC*width*a*b
    ball=4*np.pi*NC*(vel*T)**3/3
    return {'mX_t':T,'velocity':vel,'duration_mX_inverse':duration,'radius_inner':a,'radius_outer':b,
            'effective_shell_core_count':effective,'exact_poisson_cv':1/np.sqrt(effective),
            'solid_ball_core_count':ball,'incorrect_ball_cv':1/np.sqrt(ball)}

def asymptotic_check(z,dt,T):
    t,e,plus,minus=pulse_series(z,dt);edges=z['energy_edges'];de=edges[1]-edges[0];minus=minus/de
    N=np.array([simpson(F,x=t,axis=0) for F in (plus,minus)])
    total=np.array([np.trapezoid(row,x=e)+.5*de*(row[0]+row[-1]) for row in N])
    # In the T >> source-duration, T >> r0 limit, v = r/T and
    # CV^2 = integral dv [N(E(v))*gamma^3]^2 / (4 pi n_c T^3 Ntot^2).
    to_v=lambda E:np.sqrt(1-1/E**2)
    boundaries=to_v(np.r_[edges[0],e,edges[-1]])
    integrals=[]
    for row in N:
        def integrand(v):
            g=1/np.sqrt(1-v*v);return (np.interp(g,e,row)*g**3)**2
        integrals.append(sum(quad(integrand,a,b,epsabs=1e-7,epsrel=2e-9,limit=50)[0]
                             for a,b in zip(boundaries[:-1],boundaries[1:])))
    cv=np.sqrt(np.asarray(integrals)/(4*np.pi*NC*T**3))/total
    return {'mX_t':T,'density_cv':cv.tolist(),'source_total':total.tolist(),
            'scope':'Independent one-dimensional asymptotic integral for the same piecewise-linear source.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='source_t400');p.add_argument('--nv',type=int,default=4001)
    p.add_argument('--dt',type=float,default=10.);p.add_argument('--chunk',type=int,default=128)
    p.add_argument('--times',type=float,nargs='+',default=[100000.,200000.,500000.,1000000.,10000000.]);p.add_argument('--tag',required=True)
    a=p.parse_args();z=np.load(OLD/(a.source+'.npz'))
    rows=[calc(z,T,a.nv,a.dt,a.chunk) for T in a.times]
    out={'configuration':vars(a),'source_time_width':float(z['window_width_mX_t']),
         'ncore_over_mX3':NC,'rows':rows,'mono_shell_control':mono_shell(1e7,.7,200.,float(z['probe_radius_mX_inverse'])),
         'asymptotic_check':asymptotic_check(z,a.dt,1e10),
         'scope':'Synchronous identical direct-primary sources with continuous piecewise-linear energy flux, flat bin endcaps, independent homogeneous Poisson centers, and free radial streaming. Daughters, interactions, physical wave matching, actual core correlations and source variation omitted. Pointwise CV is not a collision rate or a guarantee for every point.'}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
