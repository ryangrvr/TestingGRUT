"""Frozen CP08 operational response diagnostic, no new Hamiltonian."""
import importlib.util,json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
root=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('cp08',root/'constructive_physics/attempt08/check.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
H,states,pairs,index,*_=c.build(2);H=H.toarray()
A=np.kron(np.eye(3),c.pauli(2,1,'Z').toarray())
B=np.kron(np.eye(3),c.pauli(2,0,'X').toarray())
Y=np.array([[0,-1j],[1j,0]])
comm=lambda a,b:a@b-b@a
Q=-.5j*comm(comm(H,comm(H,A)),B)[:4,:4]
assert np.linalg.norm(Q+2*np.kron(Y,Y))<1e-12
plus=np.array([1,1j])/np.sqrt(2);minus=np.array([1,-1j])/np.sqrt(2)
v=np.zeros(12,complex);v[:4]=np.kron(plus,minus)
assert abs(np.vdot(v[:4],Q@v[:4])-2)<1e-12
rows=[]
for t in [.04,.02,.01]:
 U=expm(-1j*H*t);At=U.conj().T@A@U
 R=float((1j*np.vdot(v,comm(At,B)@v)).real)
 eps=1e-5
 vals=[]
 for sign in [1,-1]:
  w=U@expm(1j*sign*eps*B)@v
  vals.append(float(np.vdot(w,A@w).real))
 kick=(vals[0]-vals[1])/(2*eps)
 assert abs(kick-R)<1e-9
 rows.append(dict(t=t,response=R,response_over_t2=R/t**2,kick_derivative=kick))
assert abs(rows[-1]['response_over_t2']-2)<.01
H0,*_=c.build(2,0);U0=expm(-1j*H0.toarray()*.4)
assert np.linalg.norm(comm(U0.conj().T@A@U0,B))<1e-12
out={'exact_projected_coefficient':'-2 Y_0 Y_1','initial_geometry':'empty order','initial_matter':'Y+ tensor Y-','controls':rows,'gamma0_cross_commutator':'zero'}
Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
print('CP08_OPERATIONAL_CROSS_RESPONSE_PASS');print(json.dumps(out,indent=2))
