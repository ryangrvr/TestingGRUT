from pathlib import Path
import sys,json
import numpy as np
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,hop
H,states,lookup,degree=build(4,4,2);mask=45;sea=sum(1<<(4+e) for e in range(6));F=[0]*12;vac=lookup[mask,sea,tuple(F)];bits,sg=hop(sea,0,4);F[0]=1;pair=lookup[mask,bits,tuple(F)]
psi=np.zeros(len(states),complex);psi[vac]=1/np.sqrt(2);psi[pair]=1j*sg/np.sqrt(2)
B=np.array([(s[1]>>0)&1 for s in states]);A=np.array([(s[1]>>2)&1 for s in states]);T=.4
w=expm_multiply(-1j*T*H,np.column_stack([psi,B*psi]));R=float(-2*np.imag(np.vdot(w[:,0],A*w[:,1])));rows=[]
for eps in [1e-3,1e-4,1e-5]:
 kick=expm_multiply(-1j*T*H,np.column_stack([psi*np.exp(1j*eps*B),psi*np.exp(-1j*eps*B)]));finite=float((A@abs(kick[:,0])**2-A@abs(kick[:,1])**2)/(2*eps));assert abs(finite-R)<1e-9;rows.append(dict(epsilon=eps,finite_impulse_response=finite,error=abs(finite-R)))
r=dict(retarded_response=R,controls=rows);Path(__file__).with_name('impulse_results.json').write_text(json.dumps(r,indent=2)+'\n');print(r)
