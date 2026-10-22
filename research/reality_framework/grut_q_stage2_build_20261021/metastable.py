#!/usr/bin/env python3
"""An explicit positive-energy version from the multiple-vacuum thread.

The small gap is supplied, with its calibration cost explicit. This is
not a mechanism selecting its magnitude or preparing the q flux.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
P=json.loads((ROOT/'INPUTS.json').read_text()); sc=P['scenario']; cfg=P['metastable_version']
lam,vev,delta=sc['lambda_s'],sc['v_s'],cfg['delta_toy_supplied']
m=cfg['m_q']; d0=sc['rho_d_initial']; r0=sc['rho_r_initial']


def potential(s):
    z=s/vev
    V=lam*(s*s-vev*vev)**2/4+delta*(2-3*z+z**3)/4
    dV=lam*s*(s*s-vev*vev)+3*delta*(z*z-1)/(4*vev)
    return V,dV


def run(s0,refined=False):
    u0,_=potential(s0)
    x0=(d0+4*u0)/(3*m*m)
    h0=math.sqrt((d0+r0+u0+.75*m*m*x0*x0)/(3*(1+x0-d0/(m*m))))
    v0=-h0*d0/(m*m)
    def rhs(t,y):
        N,x,v,h,s,w=y
        rd=d0*math.exp(-3*N); u,du=potential(s)
        return [h,v,-3*h*v-m*m*x+(rd+4*u-w*w)/3,
                .5*m*m*x-2*h*h,w,-3*h*w-du]
    ts=np.linspace(0,sc['t_end'],sc['samples'])
    sol=solve_ivp(rhs,(0,sc['t_end']),[0,x0,v0,h0,s0,0],method='DOP853',
                  t_eval=ts,rtol=1e-11 if refined else 1e-9,atol=1e-14 if refined else 1e-12,max_step=.8/m)
    if not sol.success: raise RuntimeError(sol.message)
    N,x,v,h,s,w=sol.y; u,_=potential(s)
    rd=d0*np.exp(-3*N); rr=r0*np.exp(-4*N)
    residual=3*(1+x)*h*h+3*h*v-rd-rr-.5*w*w-u-.75*m*m*x*x
    scale=3*(1+x)*h*h+np.abs(3*h*v)+rd+rr+.5*w*w+np.abs(u)+.75*m*m*x*x
    tail=ts>.8*ts[-1]
    qdec=1-.5*m*m*x/(h*h)
    return {'s_initial':s0,'s_final':float(s[-1]),'V_final':float(u[-1]),
       'H_final':float(h[-1]),'a_final':float(math.exp(N[-1])),
       'tail_average_qdec':float(np.trapezoid(qdec[tail],ts[tail])/(ts[-1]-ts[tail][0])),
       'relative_constraint_max':float(np.max(np.abs(residual)/scale))},sol


def main():
    records=[]
    for s0 in cfg['initial_s_values']:
        row,_=run(s0); records.append(row)
    _,coarse=run(cfg['initial_s_values'][1]); _,fine=run(cfg['initial_s_values'][1],True)
    refinement=float(np.max(np.abs(coarse.y[3]/fine.y[3]-1)))
    # Exact rational differentiation of the potential at its two minima.
    lv,vv,dd=Fraction(4),Fraction(1,2),Fraction(3,1000)
    mfalse2=2*lv*vv*vv-Fraction(3,2)*dd/(vv*vv)
    mtrue2=2*lv*vv*vv+Fraction(3,2)*dd/(vv*vv)
    gap_relation=(vv*vv/Fraction(3))*(mtrue2-mfalse2)-dd
    # Exact homogeneous de Sitter solution at s=-v, c_q=0.
    H=math.sqrt(delta/3); x=4*delta/(3*m*m)
    de_sitter_residual=3*(1+x)*H*H-delta-.75*m*m*x*x
    phy=cfg['physical_illustration']; numbers=P['physical_numbers']
    dp=numbers['rho_Lambda_reference_GeV4']; vp=phy['v_s_GeV_supplied']; lp=phy['lambda_s_supplied']
    mq=phy['m_q_over_M_supplied']*numbers['M_reduced_GeV']
    Hphys=math.sqrt(dp/(3*numbers['M_reduced_GeV']**2))
    physical={'v_s_GeV_supplied':vp,'lambda_s_supplied':lp,'m_q_GeV_supplied':mq,
       'delta_GeV4_calibrated':dp,'scalar_mass_at_zero_bias_GeV':math.sqrt(2*lp)*vp,
       'mass_squared_splitting_GeV2':3*dp/vp**2,
       'relative_gap_delta_over_lambda_v4':dp/(lp*vp**4),
       'false_vacuum_H_GeV':Hphys,'H_squared_over_mq_squared':(Hphys/mq)**2,
       'scope':'Analytic low-curvature illustration, not a direct integration of the physical hierarchy, collider prediction or independently calibrated benchmark.'}
    outputs={'label':cfg['label'],'gap_input_status':'Supplied, or calibrated to observed DE when physical units are attached.',
       'physical_illustration':physical,
       'runs':records,'refined_H_max_relative_difference':refinement,
       'stationary_minima':{'V_false':delta,'V_true':0.,'m_false_squared':float(mfalse2),
         'm_true_squared':float(mtrue2),'exact_gap_relation_residual':str(gap_relation)},
       'exact_false_vacuum_de_Sitter':{'H':H,'x_Fminus1':x,'constraint_residual':de_sitter_residual},
       'checks':{'positive_classical_curvatures':mfalse2>0 and mtrue2>0,
          'small_false_phase_perturbation_remains_false':records[1]['s_final']<0,
          'false_phase_accelerates':records[0]['tail_average_qdec']<0,
          'true_phase_control_decelerates':records[2]['tail_average_qdec']>0,
          'constraint_below_1e-8':all(r['relative_constraint_max']<1e-8 for r in records),
          'refinement_below_1e-6':refinement<1e-6,'exact_gap_relation':gap_relation==0},
       'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'INPUTS.json']},
       'scope':cfg['quantum_status']}
    outputs['checks']={k:bool(v) for k,v in outputs['checks'].items()}
    (ROOT/'METASTABLE.json').write_text(json.dumps(outputs,indent=2,allow_nan=False)+'\n')
    print(json.dumps(outputs,indent=2))
    if not all(outputs['checks'].values()): raise SystemExit('A metastable-version control failed.')


if __name__=='__main__':main()
