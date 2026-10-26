"""Reimplemented review case-L proxy with explicit density and step controls.

No kinetic thermalization or lump production is inferred from this proxy.
Input 4->2 and cold conversion coefficients are quoted from the review.
"""
import time
import numpy as np
from scipy.optimize import brentq
from thermo import BoseGas,M

YQ=4.8294009076988175e-9;DM=4.4e-10;MRED=2.435e18;Q=200000
UC=.8611951308217225;EC=.9233161353726057
A4=.0842146549366979;B4=.002361976125302127
SV=1.13039e-16/1.167e-17
V=.924089637375362;LP=.03;LX=.001;GX=.1
XC=np.sqrt(LP*V*V/(2*GX*GX))
RHO_C=3*LX*XC**4;NC=2*np.sqrt(2*LX)*XC**3;V0=LP*V**4/4


def gstar(T):
    ts=[1e-6,5e-4,1e-3,.1,.15,.3,1.,5.,100.,200.,1e19]
    gs=[3.36,3.36,10.75,10.75,17.25,61.75,72.25,86.25,96.25,106.75,106.75]
    return float(np.interp(np.log(T),np.log(ts),gs))


def gs(T,table='continuous'):
    if table=='review':
        x=gstar(T);return 3.91 if x<3.5 else x
    # Separate monotone entropy table removes the review interpolant's
    # 3.5 -> 3.91 jump, which can produce negative inferred time steps.
    ts=[1e-6,5e-4,1e-3,.1,.15,.3,1.,5.,100.,200.,1e19]
    ss=[3.91,3.91,10.75,10.75,17.25,61.75,72.25,86.25,96.25,106.75,106.75]
    return float(np.interp(np.log(T),np.log(ts),ss))


def entropy(T,table='continuous'):return 2*np.pi**2/45*gs(T,table)*T**3
def hubble(T):return np.sqrt(np.pi**2*gstar(T)/90)*T*T/MRED


def run(epsilon,R,ymult=1.,steps=520,nodes=96,endT=1e-8,scheme='midpoint',rate_factor=1.,case='L',table='continuous'):
    begin=time.perf_counter();gas=BoseGas(nodes);ye=YQ*ymult
    nq0=epsilon*NC;tv0=brentq(lambda T:ye*entropy(T,table)-nq0,1e-5,1e4)
    uc,ec=(UC,EC) if case=='L' else (1.,1.)
    E0=RHO_C*(1+R)+V0
    initial=gas.state(nq0,uc,energy=E0,ec=ec)
    energy_error=abs(initial['rho']+ec*M*initial['nc']-E0)/E0
    sigma=initial['s']/nq0;p=initial['nmi']/nq0
    tvs=np.geomspace(tv0,endT,steps+1)
    history=[];max_cond=0.;max_constraint=0.;last_lna=-np.log(gs(tv0,table)**(1/3)*tv0);negative_steps=0
    previous=tv0
    for i,tv in enumerate(tvs):
        lna=-np.log(gs(tv,table)**(1/3)*tv)
        te=tv if scheme=='review_endpoint' or i==0 else np.sqrt(previous*tv)
        nq=ye*entropy(te,table);state=gas.state(nq,uc,sigma=sigma)
        nc=state['nc'];peq=state['nmi']/nq;hh=hubble(te)
        max_cond=max(max_cond,nc/nq)
        max_constraint=max(max_constraint,abs(state['nq']+nc-nq)/nq,abs(state['s']/nq-sigma)/max(sigma,1.))
        def rate(pn):
            nmi=pn*nq;npl=nmi+nq-(nc if case=='L' else 0)
            gamma=rate_factor*(A4*npl*npl*nmi/4+B4*npl**3/6+B4*npl*nmi*nmi/6)
            if case=='L':gamma+=nc/Q*SV
            return gamma/hh
        if i:
            dl=lna-last_lna
            if dl<0:
                negative_steps+=1
                if table!='review':raise ValueError('entropy table moved time backward')
            if scheme=='review_endpoint':
                p=peq+(p-peq)*np.exp(-min(rate(p)*dl,700))
            elif p!=peq:
                old=p
                # Solve for a surviving fraction in [0,1]. Unlike
                # subtracting two nearly identical populations, these
                # bracket signs remain exact when the reaction is tiny.
                def equation(survival):
                    new=peq+(old-peq)*survival
                    return survival-np.exp(-min(rate((new+old)/2)*dl,700))
                survival=brentq(equation,0.,1.,xtol=2e-13,rtol=1e-11)
                p=old+(peq-old)*(1-survival)
        # A maximum possible future conversion optical depth, assuming ALL
        # net charge bound in cores and RD continuing, scales as T_vis.
        tailmax=ye*entropy(tv,table)/Q*SV/hubble(tv)
        if i%max(1,steps//130)==0 or i==steps:
            history.append(dict(Tvis=tv,Th_over_m=state['T']/M,cond_frac=nc/nq,peq=peq,p=p,
                                Gamma_over_H=rate(p),max_future_conversion_optical_depth=tailmax))
        last_lna=lna;previous=tv
    f=2*M*YQ*p/DM
    nqend=ye*entropy(endT,table);nmaxplus=(p+1)*nqend;nmaxminus=p*nqend
    future4=rate_factor*(A4*nmaxplus**2*nmaxminus/4+B4*nmaxplus**3/6+B4*nmaxplus*nmaxminus**2/6)/hubble(endT)/7
    return dict(epsilon=epsilon,R=R,ymult=ymult,case=case,steps=steps,nodes=nodes,endT=endT,scheme=scheme,
        rate_factor=rate_factor,entropy_table=table,negative_inferred_time_steps=negative_steps,
        Tvis_onset=tv0,T_hidden_initial_over_m=initial['T']/M,
        initial_cond_frac=initial['nc']/nq0,initial_energy_relative_error=energy_error,
        entropy_per_charge=sigma,p_final=p,f_extra=f,
        max_cond_frac=max_cond,final_cond_frac=nc/nq,max_thermo_constraint_error=max_constraint,
        optimistic_future_conversion_depth_bound=tailmax,
        optimistic_future_4to2_optical_depth_bound=future4,
        f_after_max_future_conversion_lower_bound=f*np.exp(-tailmax-future4),
        history=history,runtime_seconds=time.perf_counter()-begin)


if __name__=='__main__':
    import argparse,json
    a=argparse.ArgumentParser();a.add_argument('--eps',type=float,default=.15741283451037602)
    a.add_argument('--R',type=float,default=230);a.add_argument('--ymult',type=float,default=5.69009662564405)
    a.add_argument('--steps',type=int,default=520);a.add_argument('--out',required=True);s=a.parse_args()
    r=run(s.eps,s.R,s.ymult,s.steps)
    from pathlib import Path
    Path(s.out).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print({k:r[k] for k in ['epsilon','R','ymult','f_extra','max_cond_frac','entropy_per_charge','runtime_seconds']})
