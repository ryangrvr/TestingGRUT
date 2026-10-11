"""Independent JW matrix, signed permutation, propagation and asymptotic-bound controls."""
from pathlib import Path
import itertools,json,math
import numpy as np
from scipy.linalg import block_diag,expm
from check import build_joint, adjacency, energy, MASS,HOP,COST,GAMMA

# Full tensor Jordan-Wigner action, independent of sparse hop helper.
n=3;edges=list(itertools.combinations(range(n),2));M=len(edges);F=n+M
I=np.eye(2);Z=np.diag([1.,-1.]);ann=np.array([[0.,1.],[0.,0.]])
cs=[]
for j in range(F):
 op=np.array([[1.]])
 for l in reversed(range(F)):op=np.kron(op,Z if l<j else ann if l==j else I)
 cs.append(op)
H,occ,_=build_joint(n);blocks=[]
for mask in range(1<<M):
 block=MASS*M*np.eye(1<<F)
 for j,c in enumerate(cs):block+=(MASS if j<n else -MASS)*(c.T@c)
 for e,(i,j) in enumerate(edges):
  if mask>>e&1:
   block+=HOP*((cs[i].T-cs[j].T)@cs[n+e]+cs[n+e].T@(cs[i]-cs[j]))
 block=block[np.ix_(occ,occ)]+(COST*mask.bit_count()+GAMMA*M/n)*np.eye(len(occ))
 blocks.append(block)
independent=block_diag(*blocks)
for mask in range(1<<M):
 for e in range(M):
  other=mask^(1<<e);a=mask*len(occ);b=other*len(occ)
  independent[a:a+len(occ),b:b+len(occ)]-=GAMMA/n*np.eye(len(occ))
jw_error=float(np.max(abs(independent-H.toarray())));assert jw_error<1e-12

symmetry=[]
for n in [3,4]:
 H,occ,edges=build_joint(n);M=len(edges);d=len(occ);ei={e:i for i,e in enumerate(edges)};oi={b:i for i,b in enumerate(occ)}
 perm=list(range(n));perm[0],perm[1]=perm[1],perm[0]
 mp=list(perm);sg=[1]*n;edge_perm=[]
 for i,j in edges:
  a,b=perm[i],perm[j];e=ei[tuple(sorted((a,b)))];edge_perm.append(e);mp.append(n+e);sg.append(1 if a<b else -1)
 f_map=[];f_sign=[]
 for bits in occ:
  occupied=[i for i in range(n+M) if bits>>i&1];mapped=[mp[i] for i in occupied]
  inversions=sum(mapped[a]>mapped[b] for a in range(len(mapped)) for b in range(a+1,len(mapped)))
  f_map.append(oi[sum(1<<v for v in mapped)]);f_sign.append((-1)**inversions*math.prod(sg[i] for i in occupied))
 mapping=[];signs=[]
 for mask in range(1<<M):
  new=sum(1<<edge_perm[e] for e in range(M) if mask>>e&1)
  mapping.extend(new*d+np.asarray(f_map));signs.extend(f_sign)
 signs=np.asarray(signs);transformed=H[mapping,:][:,mapping].multiply(signs[:,None]).multiply(signs[None,:]).tocsr()
 err=float(np.linalg.norm((transformed-H).data));assert err<1e-12;symmetry.append(dict(N=n,error=err))

# On a path, vertex distance r requires 2r incidence hops.
n=7;B=np.zeros((n,n-1))
for e in range(n-1):B[e,e]=1;B[e+1,e]=-1
D=np.block([[np.eye(n),B],[B.T,-np.eye(n-1)]])
assert np.max(abs((D@D)[:n,:n]-(np.eye(n)+B@B.T)))<1e-12
prop=[]
for distance in [1,2,3]:
 for power in range(2*distance):assert abs(np.linalg.matrix_power(D,power)[distance,0])<1e-12
 coefficient=float(np.linalg.matrix_power(D,2*distance)[distance,0]);assert abs(coefficient)==1
 ratios=[]
 for t in [.08,.04,.02]:
  probability=abs(expm(-1j*t*D)[distance,0])**2
  leading=t**(4*distance)/math.factorial(2*distance)**2
  ratios.append(float(probability/leading))
 assert abs(ratios[-1]-1)<.001
 prop.append(dict(distance=distance,first_amplitude_power=2*distance,probability_ratios=ratios))

# Exact gap between uniform stability bound and tested full ground energies.
data=json.loads(Path(__file__).with_name('results.json').read_text())
C=((np.sqrt(2)*HOP+np.sqrt(2*HOP*HOP+2*COST*GAMMA))/(2*COST))**2
bounds=[]
for row in data['joint']:
 n=row['N'];M=n*(n-1)//2;umin=next(r['minimum'] for r in data['conditional'] if r['N']==n);a=GAMMA*M/n
 assert row['mean_edges']<=C*n
 assert umin+a-GAMMA*np.sqrt(2*C*n)<=row['energies'][0]<=umin+a+1e-10
 bounds.append(dict(N=n,static_trial=umin+a,joint_ground=row['energies'][0],mixing_gain=umin+a-row['energies'][0]))
out=dict(jordan_wigner_matrix_error=jw_error,signed_permutation=symmetry,path_propagation=prop,sparsity_constant=float(C),full_ground_bounds=bounds)
Path(__file__).with_name('adversarial_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('CP11_ADVERSARIAL_CHECKS_PASS')
