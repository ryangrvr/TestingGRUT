#!/usr/bin/env python3
"""Four stage-2 construction calculations. No old solver is imported.

See INPUTS.json and MATHEMATICS.md. Scientific scope is classical,
finite, and explicit; comparisons use already exposed reference inputs.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'INPUTS.json').read_text())
SC=P['scenario']


def potential(s):
    lam,vev=SC['lambda_s'],SC['v_s']
    return lam*(s*s-vev*vev)**2/4,lam*s*(s*s-vev*vev)


def scenario(m,c,refined=False):
    d0,r0=SC['rho_d_initial'],SC['rho_r_initial']
    s0,w0=SC['s_initial'],SC['sdot_initial']
    assert w0==0, 'Initial constraint formula below assumes zero initial scalar speed.'
    u0,_=potential(s0)
    trace0=d0+4*u0-w0*w0+4*c
    x0=trace0/(3*m*m)
    h0=math.sqrt((d0+r0+w0*w0/2+u0+c+.75*m*m*x0*x0)
                 /(3*(1+x0-d0/(m*m))))
    v0=-h0*d0/(m*m)

    def rhs(t,y):
        N,x,v,h,s,w=y
        rd=d0*math.exp(-3*N)
        u,du=potential(s)
        return [h,v,-3*h*v-m*m*x+(rd+4*u-w*w+4*c)/3,
                .5*m*m*x-2*h*h,w,-3*h*w-du]

    times=np.linspace(0,SC['t_end'],SC['samples'])
    sol=solve_ivp(rhs,(times[0],times[-1]),[0,x0,v0,h0,s0,w0],
                  method='DOP853',t_eval=times,dense_output=True,
                  rtol=1e-11 if refined else SC['rtol'],
                  atol=1e-14 if refined else SC['atol'],max_step=.8/m)
    if not sol.success: raise RuntimeError(sol.message)
    N,x,v,h,s,w=sol.y
    rd=d0*np.exp(-3*N); rr=r0*np.exp(-4*N)
    u,du=potential(s)
    rho_s=.5*w*w+u; p_s=.5*w*w-u
    constraint=3*(1+x)*h*h+3*h*v-rd-rr-rho_s-c-.75*m*m*x*x
    denominator=3*(1+x)*h*h+np.abs(3*h*v)+rd+rr+rho_s+abs(c)+.75*m*m*x*x
    trace=rd+4*u-w*w+4*c
    trace_dot=h*(-3*rd+6*w*w)+6*du*w
    xqs=trace/(3*m*m); vqs=trace_dot/(3*m*m)
    q_osc_proxy=.75*((v-vqs)**2+m*m*(x-xqs)**2)
    hdot=.5*m*m*x-2*h*h
    qdec=-1-hdot/(h*h)
    tail=times>=times[-1]*.8
    crossed=np.flatnonzero(s>=.5*SC['v_s'])
    records={
        'm_q':m,'c_supplied':c,'H_initial':h0,'N_final':float(N[-1]),
        'a_final':float(np.exp(N[-1])),'H_final':float(h[-1]),
        'scalar_first_half_vev_time':None if not len(crossed) else float(times[crossed[0]]),
        'scalar_final':float(s[-1]),'scalar_energy_initial':float(rho_s[0]),
        'scalar_energy_final':float(rho_s[-1]),'rho_r_final':float(rr[-1]),
        'rho_d_final':float(rd[-1]),'q_potential_above_c_final':float(.75*m*m*x[-1]**2),
        'q_potential_above_c_max_tail':float(np.max(.75*m*m*x[tail]**2)),
        'q_transient_energy_proxy_max':float(np.max(q_osc_proxy)),
        'F_final':float(1+x[-1]),'minimum_F':float(np.min(1+x)),
        'qdec_final_instantaneous':float(qdec[-1]),
        'qdec_tail_time_average':float(np.trapezoid(qdec[tail],times[tail])/(times[-1]-times[tail][0])),
        'scalar_w_tail_ratio_of_integrals':float(np.trapezoid(p_s[tail],times[tail])/np.trapezoid(rho_s[tail],times[tail])),
        'maximum_H_over_mq':float(np.max(h/m)),
        'maximum_relative_constraint_sum_abs':float(np.max(np.abs(constraint)/denominator)),
        'nfev':sol.nfev,
    }
    return records,sol


def gr_scenario(c,times):
    d0,r0=SC['rho_d_initial'],SC['rho_r_initial']
    def rhs(t,y):
        N,s,w=y
        u,du=potential(s)
        h=math.sqrt((d0*math.exp(-3*N)+r0*math.exp(-4*N)+.5*w*w+u+c)/3)
        return [h,w,-3*h*w-du]
    sol=solve_ivp(rhs,(times[0],times[-1]),[0,SC['s_initial'],SC['sdot_initial']],
                  method='DOP853',t_eval=times,rtol=1e-11,atol=1e-14,max_step=.05)
    if not sol.success: raise RuntimeError(sol.message)
    N,s,w=sol.y; u,_=potential(s)
    h=np.sqrt((d0*np.exp(-3*N)+r0*np.exp(-4*N)+.5*w*w+u+c)/3)
    return sol,h


def physical_bound():
    p=P['physical_numbers']; M=p['M_reduced_GeV']; rho=p['rho_Lambda_reference_GeV4']
    H=p['hbar_GeV_s']*p['H0_reference_km_s_Mpc']*1000/p['Mpc_m']
    rd=3*M*M*H*H*p['Omega_m_reference']
    rows=[]
    for me in p['unscreened_mass_controls_eV']:
        m=me*1e-9
        Uq=rd*rd/(12*M*M*m*m)
        # Full dust background correction includes K and derivative terms.
        delta_rho=3*rd*rd/(4*M*M*m*m)
        rows.append({'m_eV':me,'range_m':p['hbarc_eV_m']/me,
                     'potential_piece_over_rhoLambda':Uq/rho,
                     'full_dust_background_correction_over_rhoLambda':delta_rho/rho})
    required_U_mass=rd/math.sqrt(12*M*M*rho)
    required_H_mass=math.sqrt(3)*rd/(2*M*math.sqrt(rho))
    return {'H0_GeV':H,'rho_d_reference_GeV4':rd,'rhoLambda_over_M4':rho/M**4,
            'uniform_unscreened_controls':rows,
            'mass_if_adiabatic_potential_extrapolated_to_DE_in_H0':required_U_mass/H,
            'mass_if_adiabatic_full_background_extrapolated_to_DE_in_H0':required_H_mass/H,
            'warning':'The extrapolated masses are O(H0), outside the heavy adiabatic expansion; these are breakdown diagnostics, not fitted solutions.',
            'lab_curve_visual_range_um':[45,60],
            'lab_curve_visual_mass_meV':[p['hbarc_eV_m']/(60e-6)*1e3,p['hbarc_eV_m']/(45e-6)*1e3]}


def sm_beta(t,y):
    gy,g2,g3,yt,lam=y
    loop=16*math.pi**2
    return np.array([(41/6)*gy**3,(-19/6)*g2**3,-7*g3**3,
                     yt*(4.5*yt*yt-(17/12)*gy*gy-2.25*g2*g2-8*g3*g3),
                     24*lam*lam-6*yt**4+.375*(2*g2**4+(g2*g2+gy*gy)**2)
                     +lam*(12*yt*yt-9*g2*g2-3*gy*gy)])/loop


def multiple_vacua():
    cfg=P['vacua']; lows=np.array([cfg[k] for k in ['gY_low','g2_low','g3_low']])
    b=np.array([41/6,-19/6,-7.])
    rows=[]
    for scale in cfg['critical_scales_GeV']:
        L=math.log(scale/cfg['mu_low_GeV'])
        gstar=1/np.sqrt(1/lows**2-b*L/(8*math.pi**2))
        gy,g2,g3=gstar
        ystar=((2*g2**4+(g2*g2+gy*gy)**2)/16)**.25
        start=[gy,g2,g3,ystar,0.]
        sol=solve_ivp(sm_beta,(L,0),start,method='DOP853',rtol=1e-11,atol=1e-13)
        if not sol.success: raise RuntimeError(sol.message)
        low=sol.y[:,-1]
        rows.append({'hstar_GeV_supplied':scale,'gY_star':gy,'g2_star':g2,'yt_star_critical':ystar,
                     'lambda_star_supplied':0.,'beta_lambda_star':float(sm_beta(L,start)[-1]),
                     'yt_low_conditional_one_loop':float(low[3]),'lambda_low_conditional_one_loop':float(low[4]),
                     'gauge_roundtrip_max_abs_error':float(np.max(np.abs(low[:3]-lows))),
                     'abs_lambda_gap_if_reference_DE_from_high_vacuum_gap':4*P['physical_numbers']['rho_Lambda_reference_GeV4']/scale**4})
    # Strictly convex epsilon plus two unequal matter minima: only one is matched.
    return {'critical_boundary_flows':rows,
            'convex_counterexample':{'epsilon':'q^2/2 + C','q_from_mu':1.,'mu':1.,'C':.5,
              'matter_vacuum_energies':[0.,.01], 'total_gravitating_energies':[0.,.01]},
            'scope':'Multiple-point degeneracy is additional. No top or Higgs pole mass has been predicted; no high-order matching has been run.'}


def dissipation():
    cfg=P['dissipation']; rows=[]
    for gamma in cfg['gamma_values']:
        beta=4*gamma/9
        for vi in cfg['initial_rV']:
            def rhs(N,y):
                rv,rr,tau=y
                J=beta*rv*rv
                h=math.sqrt((rv+rr)/3)
                return [-J,-4*rr+J,1/h]
            Ns=np.linspace(0,cfg['N_end'],1001)
            sol=solve_ivp(rhs,(0,cfg['N_end']),[vi,cfg['initial_rr'],0.],
                          method='DOP853',t_eval=Ns,rtol=1e-11,atol=1e-13,max_step=.2)
            if not sol.success: raise RuntimeError(sol.message)
            rv,rr,tau=sol.y
            analytic=vi/(1+beta*vi*Ns)
            row={'gamma':gamma,'rV_initial':vi,'N_final':Ns[-1],
                 'rV_final':float(rv[-1]),'rr_final':float(rr[-1]),'tau_final':float(tau[-1]),
                 'max_relative_vacuum_flow_error':float(np.max(np.abs(rv/analytic-1))),
                 'rr_over_rV_final':float(rr[-1]/rv[-1]),
                 'bare_shifts':[]}
            for C in [-1.,0.,1.]:
                # epsilon=C+q^2+q^-2, with positive q and epsilon''>0.
                q2=((C-rv)+np.sqrt((C-rv)**2+12))/2
                q=np.sqrt(q2)
                mu=2*q-2/q**3
                q02=(C+math.sqrt(C*C+12))/2; q0=math.sqrt(q02)
                row['bare_shifts'].append({'C':C,'mu_initial':float(mu[0]),'mu_final':float(mu[-1]),
                     'mu_Minkowski_target':2*q0-2/q0**3,
                     'max_rV_reconstruction_error':float(np.max(np.abs(C-q*q+3/(q*q)-rv)))})
            rows.append(row)
    p=P['physical_numbers']; M=p['M_reduced_GeV']; robs=p['rho_Lambda_reference_GeV4']/M**4
    age_s=p['age_reference_Gyr']*1e9*p['seconds_per_year']; age_tau=age_s*M/p['hbar_GeV_s']
    natural=[]
    for gamma in cfg['gamma_values']:
        nreq=9/(4*gamma)*(1/robs-1.)
        r_asym=(2*gamma*age_tau/(3*math.sqrt(3)))**(-2/3)
        natural.append({'gamma':gamma,'efolds_to_reference_rV_from_rV1':nreq,
                        'vacuum_dominated_age_asymptote_rV':r_asym,
                        'age_asymptote_GeV4':r_asym*M**4,
                        'age_asymptote_over_reference':r_asym/robs})
    neg=cfg['negative_control']
    def neg_rhs(tau,y):
        h,rv,rr,N=y
        source=4*neg['gamma']/9*abs(h)*rv*rv
        return [-2*rr/3,-source,-4*h*rr+source,h]
    def turnaround(tau,y): return y[0]
    turnaround.terminal=True; turnaround.direction=-1
    h0=math.sqrt((neg['rV_initial']+neg['rr_initial'])/3)
    nsol=solve_ivp(neg_rhs,(0,neg['tau_limit']),[h0,neg['rV_initial'],neg['rr_initial'],0.],
                   method='DOP853',events=turnaround,rtol=1e-11,atol=1e-13,max_step=.1)
    if not nsol.success or not len(nsol.t_events[0]): raise RuntimeError('Negative-vacuum turnaround control did not finish.')
    h,rv,rr,N=nsol.y
    negative={'turnaround_tau':float(nsol.t_events[0][0]),'turnaround_N':float(N[-1]),
        'rV_initial':neg['rV_initial'],'rV_at_turnaround':float(rv[-1]),
        'maximum_constraint_absolute_residual':float(np.max(np.abs(3*h*h-rv-rr)))}
    return {'runs':rows,'negative_vacuum_control':negative,'reference_age_tau':age_tau,'finite_age_scale_diagnostic':natural,
            'positive_fixed_point':'None in the expanding positive-rV radiation-plus-vacuum branch, since d rV/dN=-(4gamma/9)rV^2.'}


def main():
    scenarios=[]; central=None; checks={}
    times=np.linspace(0,SC['t_end'],SC['samples'])
    for c in SC['c_values']:
        gr,hgr=gr_scenario(c,times)
        for m in SC['m_q_values']:
            summary,sol=scenario(m,c)
            summary['max_relative_H_vs_GR_at_same_time']=float(np.max(np.abs(sol.y[3]/hgr-1)))
            summary['N_final_minus_GR']=float(sol.y[0,-1]-gr.y[0,-1])
            scenarios.append(summary)
            if m==40 and c==.003: central=(summary,sol)
    assert central is not None
    fine,fsol=scenario(40.,.003,True)
    summary,sol=central
    refinement={'m_q':40.,'c':.003,
        'max_relative_H_difference':float(np.max(np.abs(sol.y[3]/fsol.y[3]-1))),
        'max_abs_N_difference':float(np.max(np.abs(sol.y[0]-fsol.y[0]))),
        'max_abs_scalar_difference':float(np.max(np.abs(sol.y[4]-fsol.y[4])))}
    indices=np.arange(0,len(times),SC['saved_trajectory_stride'])
    with (ROOT/'central_trajectory.csv').open('w',newline='') as file:
        writer=csv.writer(file,lineterminator='\n'); writer.writerow(['t','N','x_Fminus1','Fdot','H','s','sdot'])
        for j in indices: writer.writerow([format(float(v),'.17g') for v in [times[j],*fsol.y[:,j]]])
    bounds=physical_bound(); vacua=multiple_vacua(); relax=dissipation()
    checks={
        'sampled_Friedmann_residual_below_1e-8':all(r['maximum_relative_constraint_sum_abs']<1e-8 for r in scenarios),
        'positive_F_all_scenarios':all(r['minimum_F']>0 for r in scenarios),
        'central_refinement_H_below_1e-6':refinement['max_relative_H_difference']<1e-6,
        'central_refinement_scalar_below_1e-7':refinement['max_abs_scalar_difference']<1e-7,
        'all_classical_transitions_executed':all(r['scalar_first_half_vev_time'] is not None for r in scenarios),
        'zero_residual_tail_decelerates':all(r['qdec_tail_time_average']>0 for r in scenarios if r['c_supplied']==0),
        'positive_supplied_residual_tail_accelerates':all(r['qdec_tail_time_average']<0 for r in scenarios if r['c_supplied']>0),
        'SM_gauge_roundtrip_below_1e-8':all(r['gauge_roundtrip_max_abs_error']<1e-8 for r in vacua['critical_boundary_flows']),
        'relaxation_analytic_flow_below_1e-8':all(r['max_relative_vacuum_flow_error']<1e-8 for r in relax['runs']),
        'negative_vacuum_control_turns_around_without_cancellation':relax['negative_vacuum_control']['rV_at_turnaround']<0,
    }
    checks={k:bool(v) for k,v in checks.items()}
    output={'classification':'Builder calculation of supplied classical/phenomenological models; no empirical validation.',
       'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'INPUTS.json',Path(__file__)]},
       'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
       'scenarios':scenarios,'scenario_refinement':refinement,'bound':bounds,'vacua':vacua,'dissipation':relax,'checks':checks}
    (ROOT/'RESULTS.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':checks,'scenario_digest':[{k:r[k] for k in ['m_q','c_supplied','N_final','H_final','qdec_tail_time_average','max_relative_H_vs_GR_at_same_time']} for r in scenarios],
         'vacua':vacua['critical_boundary_flows'],'relaxation':relax['finite_age_scale_diagnostic']},indent=2))
    if not all(checks.values()): raise SystemExit('A control failed; retain the results and diagnose.')


if __name__=='__main__': main()
