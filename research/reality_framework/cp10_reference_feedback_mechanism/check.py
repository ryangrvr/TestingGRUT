"""Frozen CP10 reference-sector diagnostics; no candidate or coupling repair."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh, expm

k=Omega=.25

def lap(A): return np.diag(A.sum(axis=1))-A

def hamiltonian(A,g):
 n=len(A); m=2**n
 H=np.kron(np.eye(m),k*lap(A)).astype(float)
 for bits in range(m):
  for i in range(n):
   H[bits*n:(bits+1)*n,bits*n:(bits+1)*n]+=Omega*np.eye(n)
   other=bits^(1<<i)
   H[bits*n:(bits+1)*n,other*n:(other+1)*n]-=Omega*np.eye(n)
   H[bits*n+i,bits*n+i]+=g*(1-(1-2*((bits>>i)&1)))
 return H

rows=[]
for n in [2,3,4,5]:
 empty=np.zeros((n,n)); complete=np.ones((n,n))-np.eye(n)
 star=np.zeros((n,n));star[0,1:]=star[1:,0]=1
 bip=np.zeros((n,n));a=n//2;bip[:a,a:]=1;bip[a:,:a]=1
 graphs={'empty':empty,'star_height_two':star,'balanced_height_two':bip,'total_order':complete}
 energies={name:float(eigh(hamiltonian(A,.25),subset_by_index=[0,0])[0][0]) for name,A in graphs.items()}
 exact_empty=Omega+.25-np.sqrt(Omega**2+.25**2)
 assert abs(energies['empty']-exact_empty)<1e-12
 assert all(exact_empty-1e-12<=E<=.25+1e-12 for E in energies.values())
 assert energies['star_height_two']<=energies['total_order']+1e-12
 uniform=np.ones(n)/np.sqrt(n)
 assert all(np.linalg.norm(lap(A)@uniform)<1e-12 for A in graphs.values())
 t=.37;prop=expm(-1j*k*t*lap(complete));expected=4*np.sin(k*n*t/2)**2/n**2
 assert np.max(abs(abs(prop[0,1:])**2-expected))<1e-13
 perturb=[]
 for name,A in graphs.items():
  if name=='empty':continue # degenerate reference vacuum, separate exact solution
  coefficient=-np.trace(np.linalg.inv(k*lap(A)+2*Omega*np.eye(n)))/n
  errors=[]
  for g in [.004,.002,.001]:
   E=float(eigh(hamiltonian(A,g),subset_by_index=[0,0])[0][0])
   estimate=(E-g)/g**2;errors.append(abs(estimate-coefficient))
  assert errors[-1]<2e-5
  assert errors[-1]<errors[0]/8 # even spectrum after removing g gives O(g^2) error
  perturb.append(dict(graph=name,second_order_coefficient=float(coefficient),errors=errors))
 rows.append(dict(N=n,ground_energies=energies,empty_exact=float(exact_empty),perturbation=perturb,complete_graph_pair_probability=float(expected)))
 print(rows[-1],flush=True)
Path(__file__).with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n')
print('CP10_REFERENCE_FEEDBACK_MECHANISM_PASS')
