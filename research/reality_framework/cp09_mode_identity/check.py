from pathlib import Path
import numpy as np,json
from scipy.sparse import diags,kron,eye
from scipy.sparse.linalg import eigsh
p=Path(__file__).resolve();f=p.parents[1]/'cp09_five_event_geometry/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split('out=[]')[0],str(f),'exec'),ns)
f2=p.parents[2]/'constructive_physics/attempt08/check.py';cs={'__file__':str(f2)};exec(compile(f2.read_text().split('def test():')[0],str(f2),'exec'),cs)
rows=[]
for n in [3,4]:
 H,states,pairs=ns['build'](n);M=1<<n;E,V=eigsh(H,k=16,which='SA',tol=1e-10,v0=np.ones(H.shape[0]));idx=np.argsort(E);E=E[idx];V=V[:,idx];g=V[:,0]
 counts=np.repeat([m.bit_count()/max(1,n*(n-1)) for m in states],M)
 local=kron(eye(len(states)),cs['pauli'](n,0,'X'));collective=kron(eye(len(states)),sum(cs['pauli'](n,i,'X') for i in range(n))/n)
 probes={'relation_fraction':counts*g,'local_X0':local@g,'mean_X':collective@g}
 gap=E[1]-E[0];shell=np.abs(E-E[1])<1e-7
 weights={key:float(np.sum(abs(V[:,shell].conj().T@v)**2)) for key,v in probes.items()}
 assert weights['relation_fraction']<1e-12 and weights['mean_X']<1e-12
 assert weights['local_X0']>1e-4
 rows.append(dict(N=n,gap=float(gap),computed_gap_multiplicity=int(shell.sum()),gap_weights=weights,energies=E.tolist(),all_mode_weights={k:abs(V.conj().T@v).tolist() for k,v in probes.items()}))
 print(rows[-1]['gap_weights'],flush=True)
p.with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n');print('CP09_MODE_IDENTITY_CONTROL_PASS')
