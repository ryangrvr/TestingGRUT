"""Independent dense block-Hamiltonian test of the quadratic observable theorem;
plus exact CP12 compressed diagonal and coherent-geometry current controls."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import expm,block_diag
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,GAMMA,hop
rng=np.random.default_rng(12371);worst=0.
def herm(n):
 z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));return (z+z.conj().T)/2
for trial in range(80):
 p,q=3,5;hp,hq=herm(p),herm(q);cross=.15*(rng.normal(size=(q,p))+1j*rng.normal(size=(q,p)));b=np.linalg.norm(cross,2)
 h=np.block([[hp,cross.conj().T],[cross,hq]]);op=herm(p);oq=herm(q);op/=np.linalg.norm(op,2);oq/=np.linalg.norm(oq,2);O=block_diag(op,oq)
 for T in [.001,.1,.5,1.]:
  U=expm(-1j*T*h)[:,:p];V=expm(-1j*T*hp);delta=U.conj().T@O@U-V.conj().T@op@V
  ratio=np.linalg.norm(delta,2)/(2*b*b*T*T);assert ratio<=1+1e-6;worst=max(worst,ratio)
n=4;q=4;M=6;L=2;H,states,lookup,degree=build(n,q,L);mask=45;other=39;sea=sum(1<<(n+e) for e in range(M));fields=(0,)*(2*M)
ids=np.array([i for i,s in enumerate(states) if s[0]==mask]);block=H[ids][:,ids];bare,bs,bl,dd=build(n,q,L,masks=[mask],rewire=False)
removable=np.array([sum((s[0]>>e&1) and (s[1]>>(n+e)&1) and not(s[2][2*e] or s[2][2*e+1]) for e in range(M)) for s in bs])
expected=GAMMA*(M-q)/M*removable
assert np.max(abs((block-bare).diagonal()-expected))<1e-12
psi=np.zeros(len(states),complex);psi[lookup[mask,sea,fields]]=1/np.sqrt(2);psi[lookup[other,sea,fields]]=1j/np.sqrt(2)
flag=np.array([s[0]==mask for s in states]);current=[]
for T in [.01,.003,.001]:
 w=expm_multiply(-1j*T*H,psi);delta=float(abs(w[flag])@abs(w[flag])-.5);current.append(dict(T=T,probability_change=delta,change_over_T=delta/T))
assert abs(abs(current[-1]['change_over_T'])-GAMMA/M)<1e-4
result=dict(random_block_tests=320,max_quadratic_bound_ratio=float(worst),compressed_diagonal_identity=True,coherent_graph_current=current,exact_current_magnitude=GAMMA/M)
Path(__file__).with_name('operator_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
