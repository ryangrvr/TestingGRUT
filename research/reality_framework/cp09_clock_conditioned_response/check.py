from pathlib import Path
import numpy as np,json
from scipy.linalg import expm
from scipy.sparse import block_diag
p=Path(__file__).resolve();f=p.parents[2]/'constructive_physics/attempt09/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split("out={'cases':[]}")[0],str(f),'exec'),ns);c=ns['c']
H,states,pairs,index,*_=ns['build'](2);H=H.toarray();d=len(H);Es,V=np.linalg.eigh(H);g=V[:,0]
aset=[];bset=[]
for mask in states:
 C=np.zeros((2,2),bool)
 for k,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>k&1)
 lo=set(np.flatnonzero(~C.any(axis=0)));up=set(np.flatnonzero(~C.any(axis=1)));isolated=lo&up;lo-=isolated;up-=isolated
 zero=c.pauli(2,0,'X')*0
 bset.append(sum((c.pauli(2,int(i),'X') for i in lo),zero)/max(1,len(lo)))
 aset.append(sum((c.pauli(2,int(i),'Z') for i in up),zero)/max(1,len(up)))
A=block_diag(aset).toarray();B=block_diag(bset).toarray();eps=1e-5
# Exact finite spectral clock: supplied energies opposite to system energies.
Ec=-Es;Hc=np.diag(Ec);K=np.kron(Hc,np.eye(d))+np.kron(np.eye(d),H)
def history(initial):
 amplitudes=V.conj().T@initial
 return np.array([amplitudes[n]*V[:,n] for n in range(d)]).reshape(-1)
def conditional(psi,t):
 clock=np.exp(-1j*Ec*t)/np.sqrt(d)
 w=clock.conj()@psi.reshape(d,d)
 prob=float(np.vdot(w,w).real)
 return w/np.sqrt(prob),prob
rows=[]
hp=history(expm(1j*eps*B)@g);hm=history(expm(-1j*eps*B)@g)
assert np.linalg.norm(K@hp)<1e-12 and np.linalg.norm(K@hm)<1e-12
for tau in [.1,.4,1.]:
 wp,pp=conditional(hp,tau);wm,pm=conditional(hm,tau)
 directp=expm(-1j*H*tau)@expm(1j*eps*B)@g;directm=expm(-1j*H*tau)@expm(-1j*eps*B)@g
 assert np.linalg.norm(wp-directp)<1e-11 and np.linalg.norm(wm-directm)<1e-11
 response=float((np.vdot(wp,A@wp).real-np.vdot(wm,A@wm).real)/(2*eps))
 U=expm(-1j*H*tau);At=U.conj().T@A@U;kubo=float((1j*np.vdot(g,(At@B-B@At)@g)).real)
 assert abs(response-kubo)<1e-8
 # Change coordinate t'=2t; reciprocal clock generator and matching system parameter.
 # State conditioning at t'/2 is exactly the same event correlation.
 wc,_=conditional(hp,(2*tau)/2);assert np.linalg.norm(wc-wp)<1e-12
 rows.append(dict(clock_reading=tau,conditioned_response=response,ideal_kubo_response=kubo,state_error=float(np.linalg.norm(wp-directp)),clock_effect_probability=pp))
# Clock states at two readings are generally nonorthogonal, so not exclusive sharp times.
a=np.exp(-1j*Ec*.1)/np.sqrt(d);b=np.exp(-1j*Ec*.4)/np.sqrt(d)
out=dict(constraint_residual=float(np.linalg.norm(K@hp)),rows=rows,clock_reading_overlap_squared=float(abs(np.vdot(a,b))**2),clock_dimension=d)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_CLOCK_CONDITIONED_IDEAL_RESPONSE_PASS');print(json.dumps(out,indent=2))
