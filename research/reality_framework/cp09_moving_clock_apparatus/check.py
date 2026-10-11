from pathlib import Path
import numpy as np,json
from scipy.linalg import expm
p=Path(__file__).resolve();f=p.parents[2]/'constructive_physics/attempt09/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split("out={'cases':[]}")[0],str(f),'exec'),ns);c=ns['c']
H,states,pairs,index,*_=ns['build'](2);H=H.toarray();E,V=np.linalg.eigh(H);ground=V[:,0];d=len(H)
A=np.zeros_like(H);B=np.zeros_like(H)
for mask,k in index.items():
 if mask==1: B[k*4:(k+1)*4,k*4:(k+1)*4]=c.pauli(2,0,'X').toarray();A[k*4:(k+1)*4,k*4:(k+1)*4]=c.pauli(2,1,'Z').toarray()
 if mask==2: B[k*4:(k+1)*4,k*4:(k+1)*4]=c.pauli(2,1,'X').toarray();A[k*4:(k+1)*4,k*4:(k+1)*4]=c.pauli(2,0,'Z').toarray()
Y=np.array([[0,-1j],[1j,0]]);X=np.array([[0,1],[1,0]]);H0=np.kron(H,np.eye(2));read=np.kron(np.eye(d),X)
# Chiral clock v=1, two disjoint square interaction regions centered 0 and tau.
# Solve the autonomous PDE exactly along characteristics by piecewise exponentials.
tau=.4;eps=1e-4;zeta=1e-4
U=expm(-1j*H*tau);At=U.conj().T@A@U;ideal=float((1j*np.vdot(ground,(At@B-B@At)@ground)).real)
rows=[]
for width in [.2,.1,.05,.025]:
 vals=[]
 for sign in [1,-1]:
  state=np.kron(ground,np.array([1,0]))
  state=expm(-1j*(H0-sign*eps/width*np.kron(B,np.eye(2)))*width)@state
  state=expm(-1j*H0*(tau-width))@state
  state=expm(-1j*(H0-zeta/width*np.kron(A,Y))*width)@state
  vals.append(float(np.vdot(state,read@state).real))
 response=-(vals[0]-vals[1])/(4*eps*zeta)
 rows.append(dict(width=width,record_response=response,ideal_response=ideal,absolute_error=abs(response-ideal)))
assert rows[-1]['absolute_error']<rows[0]['absolute_error']/20
assert rows[-1]['absolute_error']<1e-4
# All clock trajectories of a packet initially wholly before both regions and
# finally wholly beyond them give identical pointer statistics: free ground phase
# before source and post-detector system unitary cancel out of pointer readout.
out=dict(ideal=ideal,source_area=eps,pointer_area=zeta,clock_velocity=1,detector_separation=tau,controls=rows)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_MOVING_CLOCK_POINTER_RESPONSE_PASS');print(json.dumps(out,indent=2))
