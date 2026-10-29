"""Finite-width s-channel exposure for time-resolved co-radial primaries.

The full Breit-Wigner includes the on-shell piece once. It is not added to the
NWA result. Energy integration is continuous over a piecewise-linear source;
the narrow pole receives adaptive quadrature instead of point sampling.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson,quad
from scipy.special import roots_legendre
from scipy.fft import rfft,irfft,next_fast_len
from local_rates import pulse_series,resonance_roots,GH,MH,GAMMA

ROOT=Path(__file__).resolve().parents[1]


def cross(Ep,Ea):
    p=np.sqrt(np.maximum(Ep*Ep-1,0));pa=np.sqrt(Ea*Ea-1)
    s=2+2*(Ep*Ea-p*pa);beta=np.sqrt(np.maximum(0,1-4/s))
    return GH**4*beta/(32*np.pi*Ep*Ea*((s-MH**2)**2+MH**2*GAMMA**2))


def energy_weights(e,Ea,nquad=8):
    """Integrate sigma*v against piecewise-linear nodal basis in E+."""
    de=np.diff(e);g,gw=roots_legendre(nquad);f=(g+1)/2
    sample=e[:-1,None]+de[:,None]*f[None,:]
    values=cross(sample,Ea);base=de[:,None]*gw[None,:]/2
    L=np.sum(base*values*(1-f)[None,:],axis=1)
    R=np.sum(base*values*f[None,:],axis=1)
    roots=[x[0] for x in resonance_roots(Ea) if e[0]<x[0]<e[-1]]
    affected=set()
    for root in roots:
        i=int(np.searchsorted(e,root)-1)
        for j in range(max(0,i-2),min(len(e)-1,i+3)):affected.add(j)
    for j in affected:
        a,b=e[j],e[j+1];knots=[a,b]
        for root in roots:
            p=np.sqrt(root*root-1);delta=MH*GAMMA*p/(2*np.sqrt(3.))
            knots.extend(x for x in root+delta*np.array([-100,-10,-1,0,1,10,100]) if a<x<b)
        knots=sorted(set(knots))
        def part(power):
            return sum(quad(lambda x:float(cross(x,Ea))*((b-x)/(b-a) if power==0 else (x-a)/(b-a)),
                            x0,x1,epsabs=1e-12,epsrel=1e-8,limit=100)[0]
                       for x0,x1 in zip(knots[:-1],knots[1:]))
        L[j],R[j]=part(0),part(1)
    weights=np.zeros(len(e));weights[:-1]+=L;weights[1:]+=R
    return weights


def run(z,dt,Rmax,nrad,amin,nquad):
    t,e,plus,minus=pulse_series(z,dt);nt=len(t);dt=t[1]-t[0];N=len(e)
    assert nt*2-1<=next_fast_len(nt*2-1)
    nn=next_fast_len(2*nt-1)
    Pf=rfft(plus,n=nn,axis=0)
    r0=float(z['probe_radius_mX_inverse']) if 'probe_radius_mX_inverse' in z else (2400+.5)*.025
    u=np.linspace(0,np.log(Rmax/r0),nrad);r=r0*np.exp(u)
    v=np.sqrt(e*e-1)/e
    total=total_nwa=omitted=0.;selected=0;details=[]
    for i,Ea in enumerate(e):
        n_a=float(simpson(minus[:,i],x=t))
        if n_a<amin:omitted+=n_a;continue
        selected+=1;va=v[i]
        mf=rfft(minus[:,i],n=nn)
        circ=irfft(Pf*np.conj(mf[:,None]),n=nn,axis=0)
        corr=np.concatenate((circ[nn-(nt-1):],circ[:nt]),axis=0)*dt
        lag0=nt-1
        shift=(r[None,:]-r0)*(1/va-1/v[:,None])/dt+lag0
        low=np.floor(shift).astype(int);fraction=shift-low
        valid=(low>=0)&(low<corr.shape[0]-1)
        idx=np.clip(low,0,corr.shape[0]-2);col=np.arange(N)[:,None]
        c=((1-fraction)*corr[idx,col]+fraction*corr[idx+1,col])*valid
        K=simpson(c/r[None,:],x=u,axis=1)/(4*np.pi*va*v)
        weights=energy_weights(e,float(Ea),nquad)
        value=float(K@weights);total+=value
        nw=0.
        for Ep,vp,coef in resonance_roots(Ea):
            if e[0]<=Ep<=e[-1]:nw+=coef*np.interp(Ep,e,K)
        total_nwa+=nw
        details.append({'Eminus':float(Ea),'Nminus_bin':n_a,'finite_s_channel_events':value,'nwa_events':nw})
    details.sort(key=lambda x:x['finite_s_channel_events'],reverse=True)
    return {'finite_width_s_channel_events':total,'nwa_events_same_interpolated_K':total_nwa,
            'finite_over_nwa':total/total_nwa if total_nwa else None,'selected_antiparticle_bins':selected,
            'omitted_antiparticle_count_below_cut':omitted,'time_step':dt,'radial_points':nrad,
            'energy_nodes':len(e),'gauss_order':nquad,'top_contributions':details[:20],
            'scope':'S-channel X Xbar -> h* -> X Xbar, Breit-Wigner counted once, direct co-radial primaries only. No h daughters, direction spread, core scatter, or collision feedback.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--dt',type=float,default=20.)
    p.add_argument('--rmax',type=float,default=46868.85);p.add_argument('--nrad',type=int,default=401)
    p.add_argument('--amin',type=float,default=1e-6);p.add_argument('--nquad',type=int,default=8);p.add_argument('--tag',required=True)
    a=p.parse_args();z=np.load(ROOT/(a.source+'.npz'))
    out={'configuration':vars(a),'result':run(z,a.dt,a.rmax,a.nrad,a.amin,a.nquad)}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
