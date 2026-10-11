from pathlib import Path
import numpy as np,json
from scipy.linalg import expm
p=Path(__file__).resolve();f=p.parents[2]/'constructive_physics/attempt09/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split("out={'cases':[]}")[0],str(f),'exec'),ns)
H,states,pairs,index,*_=ns['build'](2);H=H.toarray();energies,modes=np.linalg.eigh(H);ground=modes[:,0];d=len(H)
I=np.eye(2);X=np.array([[0,1],[1,0]],complex);Z=np.diag([1,-1]);P1=(I-Z)/2
# Clock c, system s, record r. Frozen apparatus controls omega=g=1.
F=np.diag(np.repeat([int(bool(m&(1<<pairs.index((0,1))))) for m in states],4))
Hs=np.kron(np.kron(I,H),I);Hc=np.kron(np.kron((I-X)/2,np.eye(d)),I)
Hint=np.kron(np.kron(P1,F),I-X);Ht=Hs+Hc+Hint
v=np.kron(np.kron(np.array([1,0]),ground),np.array([1,0]));E=float(np.vdot(v,Ht@v).real)
rows=[]
for t in [0.,.4,1.,2.]:
 w=expm(-1j*t*Ht)@v;assert abs(np.vdot(w,w)-1)<1e-12
 parts=[float(np.vdot(w,Q@w).real) for Q in [Hs,Hc,Hint]]
 assert abs(sum(parts)-E)<1e-11
 blocks=w.reshape(2,d,2);joint=np.sum(abs(blocks)**2,axis=1)
 # Conditional role/record agreement at a declared clock result c=1.
 clock_prob=float(joint[1].sum());agreement=0.
 for a in range(d):
  role=int(F[a,a]);agreement+=abs(blocks[1,a,role])**2
 agreement=float(agreement/clock_prob) if clock_prob>1e-12 else None
 if t>0:
  ideal=np.zeros((2,d,2),complex)
  for a in range(d):ideal[0,a,int(F[a,a])]=ground[a]
  # This comparison is a diagnostic; autonomous output need not equal ideal instant copying.
 table=np.zeros((2,2))
 for a in range(d):table[int(F[a,a]),:]+=abs(blocks[1,a,:])**2
 if clock_prob>1e-12:
  table/=clock_prob; margf=table.sum(axis=1);margr=table.sum(axis=0)
  MI=float(sum(table[a,b]*np.log(table[a,b]/(margf[a]*margr[b])) for a in range(2) for b in range(2) if table[a,b]>1e-15))
 else: MI=None
 rows.append(dict(conditional_role_record_table=table.tolist(),conditional_mutual_information_nats=MI,t=t,total_energy=sum(parts),system_energy=parts[0],clock_energy=parts[1],interaction_energy=parts[2],clock1_probability=clock_prob,record1_probability=float(joint[:,1].sum()),conditional_role_record_agreement=agreement))
# Inverse of the complete joint evolution recovers the state, records included.
w=expm(-2j*Ht)@v;inverse=np.linalg.norm(expm(2j*Ht)@w-v);assert inverse<1e-11
out=dict(controls={'N':2,'omega':1,'g':1},initial_total_energy=E,min_total_eigenvalue=float(np.linalg.eigvalsh(Ht).min()),inverse_error=float(inverse),rows=rows)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_AUTONOMOUS_APPARATUS_PASS');print(json.dumps(out,indent=2))
