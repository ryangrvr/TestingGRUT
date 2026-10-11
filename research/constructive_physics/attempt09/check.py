import importlib.util,json,math
from pathlib import Path
import numpy as np
from scipy.sparse import eye,diags,kron
from scipy.sparse.linalg import eigsh,expm_multiply
p=Path(__file__).resolve();spec=importlib.util.spec_from_file_location('cp08',p.parents[1]/'attempt08/check.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
def build(n,gamma=math.pi/4):
 H,states,pairs,index,U,edges=c.build(n,gamma)
 degree=np.zeros(len(states))
 for a,b,i,j in edges:degree[a]+=1;degree[b]+=1
 # CP08 already supplies -J transport hopping plus +Omega sum X.
 Xsum=sum(c.pauli(n,i,'X') for i in range(n));M=1<<n
 H=H+kron(diags(degree),eye(M)) + c.Omega*n*eye(H.shape[0])-2*c.Omega*kron(eye(len(states)),Xsum)
 return H.tocsr(),states,pairs,index,U,edges,degree
out={'cases':[]}
for n in [2,3,4]:
 H,states,pairs,index,U,edges,degree=build(n);M=1<<n
 vals,vec=eigsh(H,k=1,which='SA',tol=1e-10,v0=np.ones(H.shape[0]))
 assert vals[0]>-1e-9
 rng=np.random.default_rng(9);v=rng.normal(size=H.shape[0])+1j*rng.normal(size=H.shape[0]);v/=np.linalg.norm(v)
 blocks=v.reshape(-1,M)
 squared=sum(np.linalg.norm(blocks[b]-U[i,j]@blocks[a])**2 for a,b,i,j in edges)
 matter=float(np.vdot(v,c.Omega*(n*v-kron(eye(len(states)),sum(c.pauli(n,i,'X') for i in range(n)))@v)).real)
 energy=float(np.vdot(v,H@v).real);assert abs(energy-squared-matter)<1e-10
 init=np.zeros(H.shape[0],complex);init[index[0]*M]=1
 evolved=expm_multiply(-1j*H,init)
 assert abs(np.vdot(evolved,evolved)-1)<1e-10
 def totalorder(mask):return mask.bit_count()==n*(n-1)//2
 gp=c.probabilities(vec[:,0],n)
 row={'N':n,'ground_energy':float(vals[0]),'squared_difference_identity_error':abs(energy-squared-matter),'ground_geometry':c.geometry(gp,states,pairs,n),'empty_ground_probability':float(gp[index[0]]),'total_order_ground_probability':float(sum(pr for pr,m in zip(gp,states) if totalorder(m))),'t1_geometry':c.geometry(c.probabilities(evolved,n),states,pairs,n)}
 H0,*_=build(n,0);ev0=expm_multiply(-1j*H0,init);row['t1_geometry_TV_vs_gamma0']=float(sum(abs(c.probabilities(evolved,n)-c.probabilities(ev0,n)))/2)
 if n==3:
  origin=1<<pairs.index((0,2));target=sum(1<<pairs.index(e) for e in [(0,1),(1,2),(0,2)])
  iv=np.zeros_like(init);iv[index[origin]*M]=1
  controls=[]
  for t in [.02,.01,.005]:
   pr=c.probabilities(expm_multiply(-1j*t*H,iv),n)[index[target]]
   controls.append({'t':t,'endpoint_probability_over_t4':float(pr/t**4)})
  assert abs(controls[-1]['endpoint_probability_over_t4']-.75)<.002
  row['recombination_controls']=controls
 out['cases'].append(row)
# gamma=0 has exact uniform geometry x matter all-X+ zero mode for every N.
H,states,*_=build(3,0);uniform=np.ones(H.shape[0],complex)/np.sqrt(H.shape[0]);assert np.linalg.norm(H@uniform)<1e-12
out['gamma0_uniform_zero_mode_residual']=float(np.linalg.norm(H@uniform))
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_POSITIVE_CONFIGURATION_KINETIC_CHECKS_PASS');print(json.dumps(out,indent=2))
