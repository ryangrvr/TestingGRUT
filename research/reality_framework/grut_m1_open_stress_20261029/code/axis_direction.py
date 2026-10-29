"""Free spherical-Hankel outgoing/incoming diagnostic for modal traces.

Finite-time windowing and the lack of a derived quantum source mean this is
not a proof that the field has become independent particles. Exclude k<0.1,
where the in/out decomposition becomes poorly conditioned.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.special import spherical_jn,spherical_yn

ROOT=Path(__file__).resolve().parents[1]

def hankel_split(F,D,k,r,ell,direction):
    z=k*r;j=spherical_jn(ell,z);y=spherical_yn(ell,z)
    dj=k*spherical_jn(ell,z,derivative=True);dy=k*spherical_yn(ell,z,derivative=True)
    ho=j+1j*direction*y;hi=j-1j*direction*y;do=dj+1j*direction*dy;di=dj-1j*direction*dy
    den=do*hi-ho*di
    return (D*hi-F*di)/den,(F*do-D*ho)/den

def run(tag,h,lmax,kmin):
    trace=np.genfromtxt(ROOT/'.work'/(tag+'_trace.csv'),names=True,delimiter=',');t=trace['tau'];dt=t[1]-t[0];N=len(t)
    w=np.sin(np.pi*(t-t[0])/(t[-1]-t[0]));nf=4*N;freq=2*np.pi*np.fft.fftfreq(nf,dt);period=nf*dt
    radius=(int(np.floor(60/h))+.5)*h;rows=[]
    def ft(a):return dt*np.fft.fft(a*w,n=nf)
    for l in range(0,lmax+1,2):
        s=str(l);Z=ft(trace['x_l'+s]+1j*trace['y_l'+s]);DZ=ft(trace['dx_l'+s]+1j*trace['dy_l'+s]);H=ft(trace['s_minus1_l'+s]);DH=ft(trace['ds_l'+s])
        area=4*np.pi*radius**2/(.01*(2*l+1))
        for name,F,D,mass,sgn in [('plus',Z,DZ,1.,-1),('minus',Z,DZ,1.,1),('neutral',H,DH,np.sqrt(6),0)]:
            k=np.sqrt(np.maximum(freq*freq-mass*mass,0));mask=(k>kmin)&(abs(freq)<=8)
            if sgn:mask &= np.sign(freq)==sgn
            omg=freq[mask];km=k[mask];Fm=F[mask];Dm=D[mask];direction=-np.sign(omg)
            A,B=hankel_split(Fm,Dm,km,radius,l,direction)
            scale=area/period*abs(omg)/(km*radius**2)
            eo=float(np.sum(scale*abs(A)**2));ei=float(np.sum(scale*abs(B)**2))
            number=scale*abs(A)**2/abs(omg)
            angular_proxy=l*(l+1)/(km*radius)**2
            out_n=float(number.sum());num_angular=float(number@angular_proxy)
            flux=float(np.sum(area/period*(-omg)*np.imag(np.conj(Fm)*Dm)))
            rows.append({'ell':l,'species':name,'outgoing_energy_diagnostic':eo,'incoming_energy_diagnostic':ei,
                         'incoming_over_outgoing':ei/eo if eo else None,'direct_omega_gradient_flux':flux,
                         'out_minus_in_residual':eo-ei-flux,'outgoing_number_diagnostic':out_n,
                         'outgoing_number_times_l_lplus1_over_kr_squared':num_angular})
    # Check the algebra independently at known fixed k using complex amplitudes.
    errs=[]
    for l in (0,2,4,6):
        k=np.array([.4,1.2,3.]);z=k*radius;ho=spherical_jn(l,z)+1j*spherical_yn(l,z);hi=np.conj(ho)
        do=k*(spherical_jn(l,z,derivative=True)+1j*spherical_yn(l,z,derivative=True));di=np.conj(do)
        A0=1.3-.4j;B0=.2+.1j;A,B=hankel_split(A0*ho+B0*hi,A0*do+B0*di,k,radius,l,np.ones(3))
        errs.extend(abs(A-A0));errs.extend(abs(B-B0))
    proxies={}
    for species in ('plus','minus','neutral'):
        subset=[x for x in rows if x['species']==species];num=sum(x['outgoing_number_diagnostic'] for x in subset)
        proxies[species]=sum(x['outgoing_number_times_l_lplus1_over_kr_squared'] for x in subset)/num if num else None
    return {'tag':tag,'h':h,'window':'single sine taper over full recorded interval','k_min':kmin,'probe_radius':radius,
            'rows':rows,'number_weighted_angular_gradient_proxy':proxies,'synthetic_hankel_split_max_absolute_error':float(max(errs)),
            'scope':'Free-wave spherical-Hankel decomposition of finite-time field modes; source extraction is sensitive to windowing, off-shell components and propagation. It is not a two-particle angle kernel.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--tag',required=True);p.add_argument('--h',type=float,required=True);p.add_argument('--lmax',type=int,default=4);p.add_argument('--kmin',type=float,default=.1);p.add_argument('--out',required=True);a=p.parse_args()
    out=run(a.tag,a.h,a.lmax,a.kmin);(ROOT/(a.out+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
