"""Continuous-energy NWA exposure for an outward co-radial primary plume.

Unlike the first attempt, this integrates across the incoming energy bins.
Velocity-dependent crossing times and the resonant partner energy are evaluated
at Gauss-Legendre quadrature points. Source density is piecewise linear in
energy, with half-bin end plateaus so the number integral is preserved.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson
from scipy.signal import correlate,correlation_lags
from scipy.special import roots_legendre
from local_rates import pulse_series,resonance_roots

ROOT=Path(__file__).resolve().parents[1]

def run(z,dt,Rmax,nrad,nquad,cut,eminanti,threshold):
    t,e,p,m=pulse_series(z,dt);dt=float(t[1]-t[0]);edges=z['energy_edges'];de=float(np.diff(edges)[0])
    # p is already dN/(dt dE); m is dN/(dt bin). Extend half a bin at each end.
    en=np.r_[edges[0],e,edges[-1]]
    p=np.column_stack((p[:,0],p,p[:,-1]));m=np.column_stack((m[:,0],m,m[:,-1]))/de
    if threshold=='zero':
        p[:,0]=0.;m[:,0]=0.
    nt=len(t);lags=correlation_lags(nt,nt,'full')*dt
    r0=float(z['probe_radius_mX_inverse']) if 'probe_radius_mX_inverse' in z else (2400+.5)*.025
    u=np.linspace(0,np.log(Rmax/r0),nrad);r=r0*np.exp(u)
    q,qw=roots_legendre(nquad);q=(q+1)/2;qw=qw/2
    events=omitted_antiparticles=outside_roots=below_threshold=0.; root_count=0; ranking=[]
    represented=float(simpson(np.trapezoid(m,x=en,axis=1),x=t))
    for j in range(len(en)-1):
        span=en[j+1]-en[j]
        for f,w in zip(q,qw):
            Ea=en[j]+span*f
            Fa=((1-f)*m[:,j]+f*m[:,j+1])*(span*w)
            na=float(simpson(Fa,x=t))
            if Ea<eminanti:
                below_threshold+=na;continue
            if na<cut:
                omitted_antiparticles+=na;continue
            va=np.sqrt(Ea*Ea-1)/Ea
            for Ep,vp,coef in resonance_roots(Ea):
                if Ep<en[0] or Ep>en[-1]:
                    outside_roots+=na;continue
                hi=int(np.searchsorted(en,Ep));lo=max(0,hi-1);hi=min(hi,len(en)-1)
                g=0. if hi==lo else (Ep-en[lo])/(en[hi]-en[lo])
                Fp=(1-g)*p[:,lo]+g*p[:,hi]
                corr=correlate(Fp,Fa,mode='full',method='fft')*dt
                shift=(r-r0)*(1/va-1/vp)
                C=np.interp(shift,lags,corr,left=0,right=0)
                value=float(coef/(4*np.pi*va*vp)*simpson(C/r,x=u))
                events+=value;root_count+=1
                ranking.append({'Eminus':float(Ea),'Eplus_resonant':float(Ep),'quadrature_antiparticles':na,'events':value})
    ranking.sort(key=lambda row:row['events'],reverse=True)
    return {'expected_unattenuated_events':events,'represented_antiparticles':represented,
            'events_per_represented_antiparticle':events/represented,
            'omitted_antiparticle_number_below_quadrature_cut':omitted_antiparticles,
            'omitted_antiparticle_number_below_energy_test_threshold':below_threshold,
            'antiparticle_quadrature_weight_with_root_outside_plus_energy_support':outside_roots,
            'root_count':root_count,'anti_energy_gauss_order':nquad,'time_step_mX_inverse':dt,'radial_points':nrad,
            'r0_mX_inverse':r0,'rmax_mX_inverse':Rmax,'largest_contributions':ranking[:12],
            'threshold_endpoint_assumption':threshold,
            'scope':'First-order collinear primary NWA exposure, including continuous anti energy. Does not measure emitted angles or include h decay daughters, off-shell diagrams, core scattering, or collision feedback.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--dt',type=float,default=10.)
    p.add_argument('--rmax',type=float,default=46868.85);p.add_argument('--nrad',type=int,default=1201)
    p.add_argument('--nquad',type=int,default=8);p.add_argument('--cut',type=float,default=1e-9)
    p.add_argument('--eminanti',type=float,default=1.);p.add_argument('--threshold',choices=['flat','zero'],default='flat');p.add_argument('--tag',required=True)
    args=p.parse_args();z=np.load(ROOT/(args.source+'.npz'))
    result={'configuration':vars(args),'result':run(z,args.dt,args.rmax,args.nrad,args.nquad,args.cut,args.eminanti,args.threshold)}
    (ROOT/(args.tag+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
