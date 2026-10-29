"""Direct primary X/Xbar inverse-decay exposure in a co-radial ballistic plume.

The input is a positive projection of a windowed signed field-flux spectrum.
The calculation is first order in collisions and does not evolve the produced
h or its daughters. Exact collinearity is a model, not a measured angular
distribution; no interacting-cloud or final-abundance claim is made.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson,quad
from scipy.signal import correlate,correlation_lags

ROOT=Path(__file__).resolve().parents[1]
GH=.2;MH=np.sqrt(6.);GAMMA=.00018756589919939715
DELTA_Y=np.arccosh(2.)


def pulse_series(z,dt):
    centers=z['window_center_mX_t'];width=float(z['window_width_mX_t']);stop=float(centers[-1])
    times=np.arange(0,stop+dt/2,dt)
    if abs(times[-1]-stop)>1e-9:raise ValueError('time step must divide source duration')
    d=times[None,:]-centers[:,None];inside=abs(d)<width/2
    kernels=np.where(inside,np.cos(np.pi*d/width)**2,0.)
    norm=simpson(kernels,x=times,axis=1)
    kernels/=norm[:,None]
    e=z['energy_centers'];de=float(np.diff(z['energy_edges'])[0])
    plus=kernels.T@z['plus_positive']/de
    minus=kernels.T@z['minus_positive'] # bin-integrated counts/time
    return times,e,plus,minus


def resonance_roots(Ea):
    y=np.arccosh(Ea);out=[]
    for v in (y+DELTA_Y,y-DELTA_Y):
        if v>=0:
            Ep=np.cosh(v);pp=np.sinh(v)
            # Integral sigma*v over Eplus at a fixed anti E, in NWA.
            coefficient=np.pi*GH**2*pp/(4*np.sqrt(3.)*Ep*Ea)
            out.append((Ep,pp/Ep,coefficient))
    return out


def bw(Ep,Ea):
    """Direct finite-width s-channel sigma*v_Moller for collinear beams."""
    pp=np.sqrt(max(Ep*Ep-1,0));pa=np.sqrt(max(Ea*Ea-1,0))
    s=2+2*(Ep*Ea-pp*pa)
    beta=np.sqrt(max(0.,1-4/s))
    return GH**4*beta/(32*np.pi*Ep*Ea*((s-MH**2)**2+MH**2*GAMMA**2))


def finite_width_control():
    rows=[]
    for Ea in (1.15,1.4,1.8,2.5):
        # A smooth, compact positive E+ density samples the actual calculation's
        # mathematical form without depending on its windowed source.
        def density(E):return np.exp(-((E-2.8)/1.3)**2)
        roots=[r for r in resonance_roots(Ea) if 1<r[0]<8]
        split=[1.,8.]
        for E,_,_ in roots:
            # Several offsets resolve the narrow Breit-Wigner on both sides.
            p=np.sqrt(E*E-1);dsdE=2*np.sqrt(3.)/p
            w=MH*GAMMA/dsdE
            split.extend(np.clip(E+np.array([-100,-10,-1,0,1,10,100])*w,1,8))
        split=sorted(set(split))
        finite=sum(quad(lambda E:bw(E,Ea)*density(E),a,b,epsabs=1e-11,epsrel=1e-8,limit=200)[0] for a,b in zip(split[:-1],split[1:]))
        nwa=sum(coef*density(E) for E,_,coef in roots)
        rows.append({'Eminus':Ea,'finite_width':finite,'NWA':nwa,'relative_difference':finite/nwa-1})
    return rows


def run(z,timestep,Rmax,nrad,override_dt=None):
    times,e,plus,minus=pulse_series(z,timestep);dt=times[1]-times[0]
    lags=correlation_lags(len(times),len(times),'full')*dt
    r0=float(z['probe_radius_mX_inverse']) if 'probe_radius_mX_inverse' in z else (2400+.5)*.025
    u=np.linspace(0,np.log(Rmax/r0),nrad);radius=r0*np.exp(u)
    rows=[];out_total=0.;miss=0.;rank=[]
    nmax=float(simpson(minus.sum(axis=1),x=times))
    for i,Ea in enumerate(e):
        pa=np.sqrt(Ea*Ea-1);va=pa/Ea
        if va==0:continue
        n_a=float(simpson(minus[:,i],x=times))
        if n_a<=0:continue
        for Ep,vp,coef in resonance_roots(Ea):
            if Ep<e[0] or Ep>e[-1]:
                miss+=n_a
                continue
            hi=int(np.searchsorted(e,Ep));lo=max(0,hi-1);hi=min(hi,len(e)-1)
            w=0 if lo==hi else (Ep-e[lo])/(e[hi]-e[lo])
            pseries=(1-w)*plus[:,lo]+w*plus[:,hi]
            correlations=correlate(pseries,minus[:,i],mode='full',method='fft')*dt
            # Roundoff from FFT correlation can be small and negative even for
            # nonnegative inputs; preserve and report its magnitude.
            neg=float(-np.minimum(correlations,0).sum()*dt)
            shifts=(radius-r0)*(1/va-1/vp)
            c=np.interp(shifts,lags,correlations,left=0,right=0)
            value=float(coef/(4*np.pi*va*vp)*simpson(c/radius,x=u))
            out_total+=value
            rank.append({'Eminus':float(Ea),'Eplus_resonant':float(Ep),'Nminus_bin':n_a,
                         'events':value,'negative_correlation_mass':neg})
    rank.sort(key=lambda x:x['events'],reverse=True)
    return {'expected_unattenuated_events':out_total,'antiparticle_count':nmax,
            'events_per_initial_antiparticle':out_total/nmax,
            'antiparticle_bin_count_with_an_unsupported_resonant_root':miss,
            'root_count':len(rank),'largest_root_contributions':rank[:15],
            'r0_mX_inverse':r0,'rmax_mX_inverse':Rmax,'radial_points':nrad,'time_step_mX_inverse':timestep,
            'time_grid_normalization_plus':float(simpson(plus.sum(axis=1),x=times)*(e[1]-e[0])),
            'time_grid_normalization_minus':nmax,
            'scope':'Unattenuated first-order co-radial direct-primary X Xbar <-> h events. No h decay daughters, finite angular spread, core scattering, interaction feedback, or overlapping cores.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--dt',type=float,default=10.)
    p.add_argument('--rmax',type=float,default=46868.85);p.add_argument('--nrad',type=int,default=1201)
    p.add_argument('--tag',required=True);args=p.parse_args();z=np.load(ROOT/(args.source+'.npz'))
    out={'configuration':vars(args),'run':run(z,args.dt,args.rmax,args.nrad),
         'finite_width_control':finite_width_control()}
    (ROOT/(args.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
