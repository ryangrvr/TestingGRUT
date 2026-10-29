"""Finite-width primary-beam exposure with continuous incoming energy.

Correlations at adjacent antiparticle energy nodes are cached and interpolated;
their velocities and Breit-Wigner integrals are evaluated at Gauss nodes.
An optional high-order control tests incoming-energy quadrature. The positive
energy integral covers node centers only: its half-bin ends are reported.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.integrate import simpson
from scipy.special import roots_legendre
from scipy.fft import rfft,irfft,next_fast_len
from local_rates import pulse_series,resonance_roots
from fullwidth import energy_weights

ROOT=Path(__file__).resolve().parents[1]

def run(z,dt,Rmax,nrad,nquad,nplusquad,cut):
    t,e,plus,minus=pulse_series(z,dt);dt=t[1]-t[0];de=float(np.diff(z['energy_edges'])[0]);nt=len(t);N=len(e)
    ep=z['energy_edges'];ea=np.r_[ep[0],e,ep[-1]]
    # minus is per bin, convert to density and preserve half-bin plateaus.
    minus=np.column_stack((minus[:,0],minus,minus[:,-1]))/de
    nn=next_fast_len(2*nt-1);Pf=rfft(plus,n=nn,axis=0);Mf=rfft(minus,n=nn,axis=0)
    r0=float(z['probe_radius_mX_inverse']) if 'probe_radius_mX_inverse' in z else (2400+.5)*.025
    u=np.linspace(0,np.log(Rmax/r0),nrad);r=r0*np.exp(u);vp=np.sqrt(e*e-1)/e
    lag0=nt-1;col=np.arange(N)[:,None]
    q,qw=roots_legendre(nquad);q=(q+1)/2;qw=qw/2
    def makecorr(i):
        a=irfft(Pf*np.conj(Mf[:,i,None]),n=nn,axis=0)
        return np.concatenate((a[nn-(nt-1):],a[:nt]),axis=0)*dt
    C0=makecorr(0);total=nwa=omitted=0.;selected=0;start=time.monotonic()
    norms=np.trapezoid(minus,x=t,axis=0)
    for j in range(len(ea)-1):
        C1=makecorr(j+1);span=ea[j+1]-ea[j]
        if max(norms[j],norms[j+1])*span<cut:
            omitted+=.5*(norms[j]+norms[j+1])*span;C0=C1;continue
        for f,w in zip(q,qw):
            E=ea[j]+span*f;va=np.sqrt(E*E-1)/E
            shift=(r[None,:]-r0)*(1/va-1/vp[:,None])/dt+lag0
            low=np.floor(shift).astype(int);fraction=shift-low
            valid=(low>=0)&(low<C0.shape[0]-1)
            idx=np.clip(low,0,C0.shape[0]-2)
            c0=(1-fraction)*C0[idx,col]+fraction*C0[idx+1,col]
            c1=(1-fraction)*C1[idx,col]+fraction*C1[idx+1,col]
            c=((1-f)*c0+f*c1)*valid
            K=simpson(c/r[None,:],x=u,axis=1)/(4*np.pi*va*vp)
            weights=energy_weights(e,float(E),nplusquad)
            value=float(K@weights);nval=0.
            for Ep,v,coef in resonance_roots(E):
                if e[0]<=Ep<=e[-1]:nval+=coef*np.interp(Ep,e,K)
            total+=span*w*value;nwa+=span*w*nval;selected+=1
        C0=C1
        if (j+1)%50==0:print('anti intervals',j+1,'/',len(ea)-1,'sec',round(time.monotonic()-start,1),flush=True)
    # Positive endcaps missing from the E integral; no claim that their rate is bounded.
    Nt=np.trapezoid(plus,x=t,axis=0)*de
    return {'finite_width_s_channel_events':total,'nwa_events_same_interpolated_K':nwa,
            'finite_over_nwa':total/nwa if nwa else None,'antiparticle_quadrature_evaluations':selected,
            'omitted_antiparticles_below_cut':float(omitted),
            'unrepresented_positive_endcap_number_approx':float((Nt[0]+Nt[-1])/2),
            'time_step':dt,'radial_points':nrad,'energy_nodes':N,'anti_gauss_order':nquad,'plus_gauss_order':nplusquad,
            'seconds':time.monotonic()-start,
            'scope':'Unattenuated continuous-energy finite-width s-channel for point-source radial primary beams; node-center positive-energy support only, no angular spread, h daughters, core, or backreaction.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--dt',type=float,default=20.)
    p.add_argument('--rmax',type=float,default=46868.85);p.add_argument('--nrad',type=int,default=301)
    p.add_argument('--nquad',type=int,default=8);p.add_argument('--nplusquad',type=int,default=8)
    p.add_argument('--cut',type=float,default=1e-7);p.add_argument('--tag',required=True);a=p.parse_args()
    z=np.load(ROOT/(a.source+'.npz'));res=run(z,a.dt,a.rmax,a.nrad,a.nquad,a.nplusquad,a.cut)
    (ROOT/(a.tag+'.json')).write_text(json.dumps({'configuration':vars(a),'result':res},indent=2,allow_nan=False)+'\n')
    print(json.dumps(res,indent=2))

if __name__=='__main__':main()
