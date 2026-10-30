"""Direct tests of the later Gemini notes; no consciousness inference."""
from pathlib import Path
import json,numpy as np
from scipy.optimize import root
import model as m
ROOT=m.ROOT
def entropy(rho):
 e=np.linalg.eigvalsh(rho);e=e[e>1e-14];return float(-e@np.log2(e))
def reduced(rho,keep):
 dims=list(range(3));a=rho.reshape([2]*6)
 for k in sorted(set(dims)-set(keep),reverse=True):
  pos=dims.index(k);a=np.trace(a,axis1=pos,axis2=pos+len(dims));dims.remove(k)
 return a.reshape(2**len(keep),2**len(keep))
def info3(rho,classical=False):
 if classical:rho=np.diag(np.diag(rho))
 return sum(entropy(reduced(rho,[i])) for i in range(3))-sum(entropy(reduced(rho,list(ij))) for ij in [(0,1),(0,2),(1,2)])+entropy(rho)
def main():
 X=np.array([[0,1],[1,0]]);Z=np.diag([1,-1]);H=np.array([[1,1],[1,-1]])/np.sqrt(2);I=np.eye(2)
 HM=np.kron(np.kron(H,I),I);phase=np.kron(np.kron(np.diag(np.exp(-.37j*np.diag(Z))),I),I)
 def cn(target):
  out=np.zeros((8,8))
  for i in range(8):out[i^(1<<(2-target)) if i&4 else i,i]=1
  return out
 c1=cn(1);c2=cn(2);gates=[HM,c1,c2,phase,phase.conj().T,c2,c1,HM];Ws=[np.eye(8,dtype=complex)]
 for u in gates:Ws.append(u@Ws[-1])
 V=np.vstack(Ws[:8])/np.sqrt(8);C=np.zeros((64,64),complex)
 for t,u in enumerate(gates):
  a=slice(8*t,8*(t+1));j=(t+1)%8;b=slice(8*j,8*(j+1));C[a,a]+=.5*np.eye(8);C[b,b]+=.5*np.eye(8);C[b,a]-=.5*u;C[a,b]-=.5*u.conj().T
 XXX=np.kron(np.kron(X,X),X);W=Ws[4]
 def score(rho0):
  rho=W@rho0@W.conj().T;s=entropy(reduced(rho0,[1,2]));coh=float(np.trace(rho@XXX).real)
  return {'quantum_I3':info3(rho),'pointer_dephased_I3':info3(rho,True),'initial_record_entropy':s,'XXX_at_t4':coh,'J_quantum_lambda1':info3(rho)-s+coh,'J_pointer_lambda1':info3(rho,True)-s+coh}
 cases=[];inputs=[]
 for bits in range(4):
  psi=(np.eye(8)[bits]+np.eye(8)[4+(bits^3)])/np.sqrt(2);inp=W.conj().T@psi;inputs.append(inp);cases.append({'initial_record_label':f'{bits:02b}',**score(np.outer(inp,inp.conj()))})
 orig=np.eye(8)[0];original=score(np.outer(orig,orig));mixed=np.zeros((8,8));mixed[0,0]=mixed[7,7]=.5;mixed[0,7]=mixed[7,0]=.3;mixed0=W.conj().T@mixed@W
 mixedscore=score(mixed0);maxscore=score(np.eye(8)/8)
 z=m.initial();xx,yy,ss=z;fr=m.real(np.array([m.xi*xx,2*np.conj(m.LOCK)*ss*np.conj(yy),m.LOCK*yy*yy]));dr=-np.linalg.pinv(m.hessian(z),rcond=1e-12)@fr*12;dz=m.complexvec(dr)
 response=float(np.imag(dz[2]/z[2]-6*dz[0]/z[0]));eps=1e-5;delta=[]
 for R in [12*np.exp(-eps),12*np.exp(eps)]:
  u=abs(m.mu+m.LOCK*R)
  def f(a):
   x,y,s=a;return [2*m.lx*x**3-3*m.k*y*x*x+m.xi*R*x,(m.M*m.M+2*m.ly*y*y)*y-m.k*x**3-2*u*s*y,2*m.ls*(s*s-m.h*m.h/2)*s-u*y*y]
  rr=root(f,abs(z));assert rr.success;zz=np.array([rr.x[0],-rr.x[1],-rr.x[2]*np.exp(1j*np.angle(m.mu+m.LOCK*R))]);delta.append(np.angle(zz[2])-6*np.angle(zz[0]))
 fd=(delta[1]-delta[0])/(2*eps)
 checks={
  'history_isometry':np.max(abs(V.conj().T@V-np.eye(8)))<1e-13,
  'history_constraint':np.max(abs(C@V))<1e-13,
  'four_orthogonal_inputs':np.max(abs(np.array(inputs).conj()@np.array(inputs).T-np.eye(4)))<1e-13,
  'four_equal_pointer_maxima':max(abs(x['J_pointer_lambda1']-2) for x in cases)<1e-12,
  'pure_GHZ_quantum_I3_zero':max(abs(x['quantum_I3']) for x in cases)<1e-12,
  'mixed_counterexample_beats_pure_quantum_score':mixedscore['J_quantum_lambda1']>cases[0]['J_quantum_lambda1']+.3,
  'original_t4_phase_score_not_maximal':original['J_pointer_lambda1']<1.8,
  'curvature_phase_static_response_FD':abs(fd-response)<1e-8}
 p=m.P;out={'checks':{k:bool(v) for k,v in checks.items()},'passed':sum(bool(x) for x in checks.values()),'total':len(checks),'four_states':cases,'original_blank_input':original,'mixed_pointer_x06':mixedscore,'maximally_mixed':maxscore,
  'curvature_response':{'d_equilibrium_phase_d_lnR':response,'finite_difference':float(fd),'analytic_mu_over_10':m.mu/10,'scope':'Static neighboring equilibria. Demonstrates absence of a universal Planck-ratio suppression of the forward curvature response; not a computed gauge-invariant cosmological perturbation or CDI source.'},
  'h_GeV':p['heavy']['h_GeV'],'h_TeV':p['heavy']['h_GeV']/1000,
  'scope':'Direct counterexamples to uniqueness and inference claims in the later notes, using the actual pinned eight-step circuit including its t=4 phase gate. A proposed selection functional is a new boundary postulate, not evidence of a conscious universe.'}
 (ROOT/'FOLLOWUP_NOTES_CHECK.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2));assert all(checks.values())
if __name__=='__main__':main()
