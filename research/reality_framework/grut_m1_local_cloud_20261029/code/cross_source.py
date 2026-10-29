"""Independent-core contribution under an ideal synchronous Poisson model.

Campbell's second factorial moment makes particles from distinct Poisson
centers independent in the ensemble average. Rotational symmetry then makes
their relative angle uniform at any observation point. This is a model of an
ensemble average, not proof of actual mixing or an event-by-event bound.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from scipy.special import roots_legendre
from local_rates import pulse_series, GH, MH, GAMMA

ROOT=Path(__file__).resolve().parents[1]
NCORE=1.2141071422639381e-15


def kernels(energy,quad_order):
    E1=energy[:,None];E2=energy[None,:]
    p1=np.sqrt(E1*E1-1);p2=np.sqrt(E2*E2-1)
    smin=2+2*(E1*E2-p1*p2);smax=2+2*(E1*E2+p1*p2)
    nwa=np.pi*GH**2/(8*E1*E2*p1*p2)*((smin<6)&(smax>6))
    # s=6+MH*GAMMA*tan(u) absorbs the sharply peaked Lorentzian.
    width=MH*GAMMA
    umin=np.arctan((smin-6)/width);umax=np.arctan((smax-6)/width)
    nodes,weights=roots_legendre(quad_order);total=np.zeros_like(nwa)
    for x,w in zip(nodes,weights):
        u=(umax+umin)/2+(umax-umin)*x/2
        s=np.maximum(4.,6+width*np.tan(u))
        total+=w*np.sqrt(np.maximum(0.,1-4/s))
    integral=total*(umax-umin)/(2*width)
    full=GH**4*integral/(128*np.pi*E1*E2*p1*p2)
    return nwa,full


def run(z,timestep,quad_order,end_times):
    t,e,plus,minus=pulse_series(z,timestep)
    de=float(np.diff(z['energy_edges'])[0]);plus*=de
    Cplus=cumulative_trapezoid(plus,x=t,axis=0,initial=0.)
    Cminus=cumulative_trapezoid(minus,x=t,axis=0,initial=0.)
    nwa,full=kernels(e,quad_order);rows=[]
    for label,K in [('narrow_width',nwa),('finite_width_s_channel',full)]:
        rate=NCORE*np.einsum('ti,ij,tj->t',Cplus,K,Cminus,optimize=True)
        early=float(simpson(rate,x=t));last=float(rate[-1]);vals=[]
        for T in end_times:
            if T<t[-1]:raise ValueError('end time must follow the final source window')
            vals.append({'mX_t':T,'events_per_reference_core':early+(T-t[-1])*last,
                         'optical_time_ncore_times_mX_t':NCORE*T})
        rows.append({'kernel':label,'during_source_events':early,'instantaneous_final_rate_per_mX_t':last,
                     'upper_bound_rate_times_time':[{ 'mX_t':T,'events_per_core':last*T} for T in end_times],
                     'integrated':vals})
    return {'source_duration_mX_t':float(t[-1]),'source_count_plus':float(Cplus[-1].sum()),
            'source_count_minus':float(Cminus[-1].sum()),'energy_bins':len(e),'rows':rows}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='source_t400');p.add_argument('--dt',type=float,default=10.)
    p.add_argument('--quad',type=int,default=64);p.add_argument('--tag',required=True);a=p.parse_args()
    result=run(np.load(ROOT/(a.source+'.npz')),a.dt,a.quad,[10000.,100000.,200000.,1000000.])
    out={'configuration':vars(a),'ncore_over_mX3':NCORE,'result':result,
         'scope':'Ensemble-mean first-order X Xbar s-channel events from distinct centers only, with identical synchronous direct-primary radial sources, homogeneous Poisson centers, vacuum straight rays, no post-10000 source, no h decay daughters, core interactions, or cosmological expansion. K includes its pole once. Not an event-by-event bound, validated mixing, or a loss count.'}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
