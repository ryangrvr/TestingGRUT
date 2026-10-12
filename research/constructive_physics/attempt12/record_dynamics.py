"""Gauge-invariant local preparation, connection locking and reversible dynamics."""
from pathlib import Path
import json,itertools
import numpy as np
from scipy.sparse.linalg import expm_multiply
from model import build,hop,ELECTRIC,MASS,HOP,GAMMA,COST

def run(L):
 n=4;q=3;M=6;mask=41;e=0;vertex=0
 H,states,lookup,degree=build(n,q,L)
 bits=sum(1<<(n+k) for k in range(M));fields=(0,)*(2*M);vac=lookup[mask,bits,fields]
 newbits,sign=hop(bits,vertex,n+e);F=list(fields);F[0]=1;pair=lookup[mask,newbits,tuple(F)]
 removed=np.array([not(s[0]>>e&1) for s in states]);changed=np.array([s[0]!=mask for s in states])
 out=[]
 for angle in [0.,np.pi/4,np.pi/2]:
  initial=np.zeros(len(states),complex);initial[vac]=np.cos(angle);initial[pair]=1j*sign*np.sin(angle)
  energy=float(np.vdot(initial,H@initial).real)
  rows=[]
  for t in [.04,.02,.01,.005]:
   w=expm_multiply(-1j*t*H,initial);p=abs(w)**2
   assert abs(np.linalg.norm(w)-1)<1e-11;assert abs(np.vdot(w,H@w).real-energy)<1e-10
   lead=(M-q)*(GAMMA/M)**2*np.cos(angle)**2
   if abs(np.cos(angle))<1e-12:
    predicted=(M-q)*(GAMMA/M)**2*HOP**2/4;ratio=float(p[removed].sum()/t**4/predicted)
   else:ratio=float(p[removed].sum()/t**2/lead)
   rows.append(dict(t=t,removed_probability=float(p[removed].sum()),geometry_departure=float(p[changed].sum()),ratio_to_leading=ratio))
  assert abs(rows[-1]['ratio_to_leading']-1)<.015
  back=expm_multiply(1j*.1*H,expm_multiply(-1j*.1*H,initial));assert np.linalg.norm(back-initial)<1e-10
  out.append(dict(angle=angle,initial_total_energy=energy,controls=rows))
 result=dict(cutoff=L,dimension=len(states),record_bare_energy=2*MASS+ELECTRIC/2,rows=out)
 Path(__file__).with_name(f'record_dynamics_L{L}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));return result
if __name__=='__main__':
 a=run(2);b=run(3)
 for x,y in zip(a['rows'],b['rows']):
  for s,t in zip(x['controls'],y['controls']):assert abs(s['removed_probability']-t['removed_probability'])<1e-12
 print('CP12_LOCAL_RECORD_DYNAMICS_PASS')
