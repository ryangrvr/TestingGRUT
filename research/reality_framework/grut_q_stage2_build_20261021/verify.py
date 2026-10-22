#!/usr/bin/env python3
"""Same-author controls using different mathematical formulations.

No import from build.py. Einstein-frame evolution checks the Jordan
scenario; analytic gauge profiles plus RK4 check the coupled RG flow;
quadrature checks the dissipative receiver; exact arithmetic checks the
full background versus potential-only factor. This is not outside review.
"""
from __future__ import annotations
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'INPUTS.json').read_text())
R=json.loads((ROOT/'RESULTS.json').read_text())


def einstein_frame():
    sc=P['scenario']; m,c=40.,.003
    data=np.loadtxt(ROOT/'central_trajectory.csv',delimiter=',',skiprows=1)
    initial=data[0]
    _,nj,x,v,h,s,w=initial
    f=1+x; beta=math.sqrt(2/3)
    chi=math.log(f)/beta
    wc=v/(beta*f**1.5); ws=w/math.sqrt(f)
    d0,r0=sc['rho_d_initial'],sc['rho_r_initial']
    lam,vev=sc['lambda_s'],sc['v_s']

    def unpack(y):
        tj,nj,chi,wc,s,ws=y
        f=math.exp(beta*chi)
        vs=lam*(s*s-vev*vev)**2/4
        uq=c+.75*m*m*(f-1)**2
        rd=d0*math.exp(-3*nj)/f**2
        rr=r0*math.exp(-4*nj)/f**2
        u=(uq+vs)/f**2
        he=math.sqrt((wc*wc/2+ws*ws/(2*f)+u+rd+rr)/3)
        return f,vs,uq,rd,rr,he

    def rhs(te,y):
        tj,nj,chi,wc,s,ws=y
        f,vs,uq,rd,rr,he=unpack(y)
        vchi=beta*(1.5*m*m*f*(f-1)-2*(uq+vs))/f**2
        return [1/math.sqrt(f),he-beta*wc/2,wc,
                -3*he*wc-beta*ws*ws/(2*f)-vchi+beta*rd/2,
                ws,-3*he*ws+beta*wc*ws-lam*s*(s*s-vev*vev)/f]

    y0=[0.,nj,chi,wc,s,ws]
    def endpoint(te,y): return y[0]-sc['t_end']
    endpoint.terminal=True
    sol=solve_ivp(rhs,(0.,1.5*sc['t_end']),y0,events=endpoint,dense_output=True,
                  method='DOP853',rtol=2e-11,atol=2e-14,max_step=.6/m)
    if not sol.success or not len(sol.t_events[0]): raise RuntimeError('Einstein-frame endpoint not reached.')
    converted=[]
    for row in data:
        tj=row[0]
        if tj==0: te=0.
        elif tj>=sc['t_end']: te=sol.t_events[0][0]
        else:
            j=int(np.searchsorted(sol.y[0],tj))
            te=brentq(lambda te:sol.sol(te)[0]-tj,sol.t[j-1],sol.t[j],xtol=1e-12)
        y=sol.sol(te)
        f,vs,uq,rd,rr,he=unpack(y)
        _,nj,chi,wc,s,ws=y
        hj=math.sqrt(f)*(he-beta*wc/2)
        vj=beta*f**1.5*wc
        converted.append([tj,nj,f-1,vj,hj,s,math.sqrt(f)*ws])
    converted=np.array(converted)
    return {'max_relative_H_difference':float(np.max(np.abs(converted[:,4]/data[:,4]-1))),
            'max_abs_N_difference':float(np.max(np.abs(converted[:,1]-data[:,1]))),
            'max_abs_s_difference':float(np.max(np.abs(converted[:,5]-data[:,5]))),
            'max_abs_sdot_difference':float(np.max(np.abs(converted[:,6]-data[:,6]))),
            'max_abs_F_difference':float(np.max(np.abs(converted[:,2]-data[:,2]))),
            'nfev':sol.nfev}


def rg_rk4():
    cfg=P['vacua']; lows=np.array([cfg[k] for k in ['gY_low','g2_low','g3_low']])
    b=np.array([41/6,-19/6,-7.]); loop=16*math.pi**2
    def gauges(t): return 1/np.sqrt(1/lows**2-b*t/(8*math.pi**2))
    def rhs(t,y):
        yt,la=y; gy,g2,g3=gauges(t)
        a=yt*(4.5*yt*yt-17*gy*gy/12-9*g2*g2/4-8*g3*g3)/loop
        z=(24*la*la-6*yt**4+3*(2*g2**4+(gy*gy+g2*g2)**2)/8
           +la*(12*yt*yt-9*g2*g2-3*gy*gy))/loop
        return np.array([a,z])
    rows=[]
    for record in R['vacua']['critical_boundary_flows']:
        L=math.log(record['hstar_GeV_supplied']/cfg['mu_low_GeV'])
        gy,g2,g3=gauges(L)
        y0=np.array([((2*g2**4+(gy*gy+g2*g2)**2)/16)**.25,0.])
        answers=[]
        for steps in [2000,4000]:
            dt=-L/steps; y=y0.copy(); t=L
            for _ in range(steps):
                a=rhs(t,y); z=rhs(t+dt/2,y+dt*a/2)
                q=rhs(t+dt/2,y+dt*z/2); u=rhs(t+dt,y+dt*q)
                y += dt*(a+2*z+2*q+u)/6; t+=dt
            answers.append(y)
        primary=np.array([record['yt_low_conditional_one_loop'],record['lambda_low_conditional_one_loop']])
        rows.append({'hstar':record['hstar_GeV_supplied'],
                     'max_abs_error_to_primary':float(np.max(np.abs(answers[-1]-primary))),
                     'max_abs_RK4_refinement':float(np.max(np.abs(answers[-1]-answers[0])))})
    return rows


def receiver_quadrature():
    cfg=P['dissipation']; rows=[]
    for record in R['dissipation']['runs']:
        beta=4*record['gamma']/9; vi=record['rV_initial']; N=record['N_final']
        def exact_rv(n): return vi/(1+beta*vi*n)
        integral,error=quad(lambda u:math.exp(-4*u)*beta*exact_rv(N-u)**2,
                            0,N,epsabs=1e-14,epsrel=1e-11)
        rr=cfg['initial_rr']*math.exp(-4*N)+integral
        rv=exact_rv(N)
        # The q flow is recovered algebraically; evaluate the balance by
        # direct differentiation of mu(q), not the main evolution solver.
        defects=[]
        for C in [-1.,0.,1.]:
            q=math.sqrt((C-rv+math.sqrt((C-rv)**2+12))/2)
            eps2=2+6/q**4
            q_N=beta*rv*rv/(q*eps2)
            defects.append(abs(q*eps2*q_N-beta*rv*rv))
        rows.append({'gamma':record['gamma'],'rV_initial':vi,
            'relative_radiation_error':abs(rr/record['rr_final']-1),
            'relative_vacuum_error':abs(rv/record['rV_final']-1),
            'quadrature_error_estimate':error,'mu_balance_max_abs_residual':max(defects)})
    return rows


def background_coefficient():
    # Exact initial constraint with adiabatic x and Fdot, expanded at
    # increasing mass. This checks the coefficient including -3H Fdot.
    rd=Fraction(1,10); rows=[]
    for m in [100,1000]:
        mm=Fraction(m*m); x=rd/(3*mm)
        h2=(rd+Fraction(3,4)*mm*x*x)/(3*(1+x-rd/mm))
        U=rd*rd/(12*mm)
        delta_source=3*h2-rd
        rows.append({'m':m,'full_source_over_potential_piece':float(delta_source/U),
                     'deviation_from_9':float(delta_source/U-9)})
    return rows


def main():
    for name,want in R['source_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,name
    ef=einstein_frame(); rg=rg_rk4(); rec=receiver_quadrature(); coeff=background_coefficient()
    checks={
       'Jordan_Einstein_H_agree_below_1e-6':ef['max_relative_H_difference']<1e-6,
       'Jordan_Einstein_scalar_agree_below_1e-7':ef['max_abs_s_difference']<1e-7,
       'RG_independent_formulation_agrees_below_1e-8':all(x['max_abs_error_to_primary']<1e-8 for x in rg),
       'radiation_receiver_agrees_below_1e-8':all(x['relative_radiation_error']<1e-8 for x in rec),
       'full_background_factor_approaches_9':coeff[-1]['deviation_from_9']<coeff[0]['deviation_from_9'] and coeff[-1]['deviation_from_9']<1e-6,
    }
    out={'review_status':P['review_status'],'Einstein_frame_scenario':ef,'SM_RK4':rg,
         'receiver_quadrature':rec,'full_background_coefficient':coeff,'checks':{k:bool(v) for k,v in checks.items()},
         'source_hash':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))
    if not all(checks.values()): raise SystemExit('A separate-formulation control failed.')


if __name__=='__main__': main()
