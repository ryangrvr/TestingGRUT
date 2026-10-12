"""Broad classical uncertainty versus genuine initial geometry coherence."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,GAMMA
H,states,lookup,degree=build(4,4,2);sea=sum(1<<(4+e) for e in range(6));E=(0,)*12;G,F=45,39
v=np.zeros((len(states),3),complex);v[lookup[G,sea,E],0]=1;v[lookup[F,sea,E],1]=1;v[:,2]=(v[:,0]+1j*v[:,1])/np.sqrt(2)
flag=np.array([s[0]==G for s in states]);b2=GAMMA**2*4*2/36;rows=[]
for T in [.01,.003,.001]:
 w=expm_multiply(-1j*T*H,v);p=(abs(w[flag])**2).sum(axis=0);mix=float((p[0]+p[1])/2);coh=float(p[2])
 assert abs(mix-.5)<=2*b2*T*T+1e-12
 rows.append(dict(T=T,mixed_change=mix-.5,coherent_change=coh-.5,mixed_quadratic_bound=2*b2*T*T,coherent_linear_coefficient=(coh-.5)/T))
neutral_ids=[i for i,s in enumerate(states) if s[1]==sea and not any(s[2])]
uniform=np.zeros(len(states));uniform[neutral_ids]=1/np.sqrt(len(neutral_ids))
coherent_energy=float(uniform@(H@uniform));dephased_energy=float(np.mean(H.diagonal()[neutral_ids]));energy_cost=GAMMA*4*2/6
assert abs(dephased_energy-coherent_energy-energy_cost)<1e-12
result=dict(rows=rows,neutral_dephasing_energy_cost=dephased_energy-coherent_energy,exact_energy_cost=energy_cost,coherence_trace_norm=1.,qualification='Dephased mixture is a supplied preparation, not dynamically derived decoherence.')
Path(__file__).with_name('mixture_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
