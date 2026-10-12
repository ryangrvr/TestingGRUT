"""Actual connectivity creation from triangle+isolate; no claim of spatial selection."""
from pathlib import Path
import sys,json,itertools
import numpy as np
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,hop,GAMMA
from thermal_control import graph_class
rows=[]
for L in [4,5]:
 n=4;q=3;M=6;mask=11;H,states,lookup,degree=build(n,q,L);sea=sum(1<<(n+e) for e in range(M));F=[0]*12;vac=lookup[mask,sea,tuple(F)];bits,sg=hop(sea,0,4);F[0]=1;pair=lookup[mask,bits,tuple(F)]
 psi=np.zeros(len(states),complex);psi[vac]=1/np.sqrt(2);psi[pair]=1j*sg/np.sqrt(2)
 B=np.array([(s[1]>>0)&1 for s in states]);A=np.array([(s[1]>>3)&1 for s in states]);conn=np.array([graph_class(n,s[0])[0] for s in states]);same=np.array([s[0]==mask for s in states]);boundary=np.array([max(abs(e) for e in s[2])==L for s in states]);w=np.column_stack([psi,B*psi]);last=0;values=[];energy=float(np.vdot(psi,H@psi).real);b2=GAMMA**2*q*(M-q)/M**2
 for T in [.4,1.,2.,4.,20.,40.]:
  w=expm_multiply(-1j*(T-last)*H,w);last=T;p=abs(w[:,0])**2;connected=float(p@conn);depart=float(1-p[same].sum());R=float(-2*np.imag(np.vdot(w[:,0],A*w[:,1])))
  assert connected<=depart+1e-10 and depart<=b2*T*T+1e-9
  assert abs(np.vdot(w[:,0],H@w[:,0]).real-energy)<1e-8
  assert abs(np.linalg.norm(w[:,0])-1)<1e-9
  values.append(dict(T=T,connected_probability=connected,departure=depart,isolated_detector_response=R,cutoff_boundary=float(p@boundary)))
 rows.append(dict(N=n,Q=q,L=L,initial_graph='triangle plus isolated vertex',energy=energy,rows=values))
for a,b in zip(rows[0]['rows'],rows[1]['rows']):assert abs(a['connected_probability']-b['connected_probability'])<1e-7
Path(__file__).with_name('assembly_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
