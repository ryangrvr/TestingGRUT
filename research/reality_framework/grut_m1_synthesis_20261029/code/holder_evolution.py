"""Evolve a candidate phase holder in a prescribed reheating background.

The homogeneous inflaton and radiation are evolved self-consistently, but the
hidden fields are spectators (their backreaction and gravitational constraint
corrections are not solved). This is a controlled comparison, not an FRW or
quantum completion of M1. Curvature is R=2 beta^2 I^2 - Idot^2 and can change
sign during inflaton oscillation.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
import phase_holder as ph
ROOT=Path(__file__).resolve().parents[1]
BETA=3.;GAMMA=.5

def unpack(v):
 I,J,rad,lna=v[:4];z=v[4:10:2]+1j*v[5:10:2];vel=v[10:16:2]+1j*v[11:16:2]
 return I,J,rad,lna,z,vel

def setup(z):
 v=np.zeros(17);v[0]=np.sqrt(6)/BETA;v[4:10:2]=z.real;v[5:10:2]=z.imag;return v

def rhs(t,v,xi):
 I,J,rad,lna,z,vel=unpack(v);H=np.sqrt((J*J/2+BETA*BETA*I*I/2+rad)/3);R=2*BETA**2*I**2-J**2
 Rdot=6*BETA**2*I*J+2*(3*H+GAMMA)*J**2
 V,F,_,_=ph.pot_force(z,xi,R);acc=-3*H*vel-F
 out=np.zeros_like(v);out[0]=J;out[1]=-(3*H+GAMMA)*J-BETA**2*I;out[2]=GAMMA*J*J-4*H*rad;out[3]=H
 out[4:10:2]=vel.real;out[5:10:2]=vel.imag;out[10:16:2]=acc.real;out[11:16:2]=acc.imag
 dVdR=ph.xiX*abs(z[0])**2+2*xi/ph.ell**5*np.real(1j*np.conj(z[2])*z[0]**6)
 out[16]=-6*H*np.vdot(vel,vel).real+Rdot*dVdR
 return out

def observables(t,V,xi):
 I,J,rad,lna,z,vel=unpack(V);H=np.sqrt((J*J/2+BETA*BETA*I*I/2+rad)/3);R=2*BETA**2*I**2-J**2
 charges=-2*np.exp(3*lna)*np.imag(np.conj(z)*vel);E=np.vdot(vel,vel).real+ph.pot_force(z,xi,R)[0]
 return [t,np.exp(lna),H,R,charges[0],charges[0]+3*charges[1],charges@[1,3,6],E,E-V[16],abs(z[0]),abs(z[1]),abs(z[2]),ph.wrap(np.angle(z[2])-6*np.angle(z[0]))]

def run(name,z,xi,end,rtol):
 v0=setup(z);ts=np.linspace(0,end,int(end*50)+1)
 def runaway(t,y):return 300-max(abs(y[4:10]));
 runaway.terminal=True;runaway.direction=-1
 tic=time.monotonic();sol=solve_ivp(lambda t,v:rhs(t,v,xi),(0,end),v0,method='DOP853',rtol=rtol,atol=rtol*1e-3,max_step=.08/BETA,t_eval=ts,events=runaway)
 if not sol.success or sol.t[-1]<end-1e-8:raise RuntimeError(f'{name} incomplete at {sol.t[-1]} {sol.message}')
 a=np.array([observables(t,v,xi) for t,v in zip(sol.t,sol.y.T)]);g=a[:,6];E0=a[0,7];ledger=a[:,8]
 mask=a[:,0]>=.8*end;mean=float(np.trapezoid(a[mask,4],a[mask,0])/(a[mask,0][-1]-a[mask,0][0]));mask2=a[:,0]>=.9*end
 mean2=float(np.trapezoid(a[mask2,4],a[mask2,0])/(a[mask2,0][-1]-a[mask2,0][0]))
 keep=np.linspace(0,len(a)-1,min(len(a),4001)).round().astype(int)
 return {'name':name,'xi':xi,'end_t_Hi':end,'initial_fields_over_H':[[float(v.real),float(v.imag)] for v in z],
   'late_mean_NX':mean,'late_last_half_mean_NX':mean2,'relative_late_window_change':abs(mean2-mean)/max(abs(mean),1e-20),
   'sampled_NX_at_20':float(a[np.argmin(abs(a[:,0]-20)),4]),'negative_ricci_sample_fraction':float(np.mean(a[:,3]<0)),
   'max_abs_gauss_charge':float(np.max(abs(g))),'energy_ledger_rel_drift':float(np.max(abs(ledger-ledger[0]))/max(abs(E0),1)),
   'initial_hidden_energy_over_H4':float(E0),'initial_hidden_over_total_background':float(abs(E0)*(ph.H/2.435e18)**2/3),
   'nfev':sol.nfev,'seconds':time.monotonic()-tic},a[keep]

def main():
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--end',type=float,default=200);p.add_argument('--rtol',type=float,default=1e-9)
 p.add_argument('--case',default='all');p.add_argument('--out',default='holder_evolution');a=p.parse_args()
 phase=json.loads((ROOT/'phase_holder.json').read_text());states=phase['states'];find=lambda xi:np.array([complex(*x) for x in next(v for v in states if v['xi']==xi)['fields_over_H']])
 cases=[('original_prepared_no_holder',ph.z0,0.),('aligned_equilibrium_no_holder',find(0.),0.),('local_equilibrium_holder_2e8',find(2e8),2e8),('local_equilibrium_holder_4e8',find(4e8),4e8)]
 if a.case!='all':
  cases=[x for x in cases if x[0]==a.case]
  if not cases:raise ValueError('Unknown case '+a.case)
 summaries=[];arrays={}
 for name,z,xi in cases:
  print('start',name,flush=True);s,arr=run(name,z,xi,a.end,a.rtol);summaries.append(s);arrays[name]=arr;print('done',name,s['seconds'],s['late_mean_NX'],flush=True)
 out={'columns':['t_Hi','a','H_over_Hi','R_over_Hi2','NX','Nlight','Gauss','hidden_energy','hidden_energy_minus_work','absX','absY','absSigma','relative_phase'],
      'background':{'beta':BETA,'Gamma_over_Hi':GAMMA,'spectator_approximation':True},'cases':summaries,
      'scope':'Prescribed homogeneous inflaton/radiation plus spectator hidden scalar and curvature phase operator. Local de Sitter equilibrium is not a state-selection law, UV completion, full metric response, quantum perturbation calculation or final abundance.'}
 np.savez_compressed(ROOT/(a.out+'.npz'),**arrays);(ROOT/(a.out+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
