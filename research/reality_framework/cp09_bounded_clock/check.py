from pathlib import Path
import numpy as np,json
from scipy.linalg import expm
p=Path(__file__).resolve();f=p.parents[2]/'constructive_physics/attempt09/check.py';ns={'__file__':str(f)};exec(compile(f.read_text().split("out={'cases':[]}")[0],str(f),'exec'),ns);c=ns['c'];H,states,pairs,*_=ns['build'](2);H=H.toarray();E,V=np.linalg.eigh(H);d=len(H);A=np.zeros_like(H);B=np.zeros_like(H)
for k,mask in enumerate(states):
 if mask in [1,2]:
  i,j=(0,1) if mask==1 else (1,0);B[k*4:(k+1)*4,k*4:(k+1)*4]=c.pauli(2,i,'X').toarray();A[k*4:(k+1)*4,k*4:(k+1)*4]=c.pauli(2,j,'Z').toarray()
A=V.conj().T@A@V;B=V.conj().T@B@V;Y=np.array([[0,-1j],[1j,0]]);X=np.array([[0,1],[1,0]]);Ak=np.kron(A,Y);Bk=np.kron(B,np.eye(2));Ek=np.repeat(E,2)
tau=.4;width=.15;eps=zeta=.001;T=1.5
U=np.diag(np.exp(-1j*E*tau));g=np.eye(1,d,0).ravel();At=U.conj().T@A@U;ideal=float((1j*np.vdot(g,(At@B-B@At)@g)).real)
def run(mass,grid,dt,sign):
 L=4.;dx=L/grid;q=np.arange(grid)*dx-L/2;mom=2*np.pi*np.fft.fftfreq(grid,d=dx)
 packet=np.exp(-(q+.7)**2/(4*.08**2)+1j*mass*q);packet/=np.linalg.norm(packet)
 state=np.zeros((grid,2*d),complex);state[:,0]=packet
 source=abs(q)<width/2;det=abs(q-tau)<width/2
 S=expm(1j*sign*eps/width*Bk*dt);D=expm(1j*zeta/width*Ak*dt)
 free=np.exp(-1j*dt/2*(mom[:,None]**2/(2*mass)+Ek[None,:]))
 def energy(w):
  pw=np.fft.fft(w,axis=0,norm='ortho');val=np.sum(abs(pw)**2*(mom[:,None]**2/(2*mass)+Ek[None,:]))
  val-=sign*eps/width*np.einsum('qi,ij,qj->',w[source].conj(),Bk,w[source]).real
  val-=zeta/width*np.einsum('qi,ij,qj->',w[det].conj(),Ak,w[det]).real
  return float(val)
 ein=energy(state)
 for _ in range(round(T/dt)):
  state=np.fft.ifft(np.fft.fft(state,axis=0,norm='ortho')*free,axis=0,norm='ortho')
  state[source]=state[source]@S.T;state[det]=state[det]@D.T
  state=np.fft.ifft(np.fft.fft(state,axis=0,norm='ortho')*free,axis=0,norm='ortho')
 assert abs(np.linalg.norm(state)-1)<1e-10
 read=float(np.einsum('qi,ij,qj->',state.conj(),np.kron(np.eye(d),X),state).real)
 return read,energy(state)-ein,ein
rows=[]
for mass,grid,dt in [(20,256,.005),(80,512,.005),(80,512,.0025),(80,1024,.0025),(320,1024,.0025),(320,2048,.00125)]:
 a=run(mass,grid,dt,1);b=run(mass,grid,dt,-1);response=-(a[0]-b[0])/(4*eps*zeta)
 rows.append(dict(mass=mass,grid=grid,dt=dt,response=response,ideal=ideal,error=abs(response-ideal),max_energy_drift=max(abs(a[1]),abs(b[1])),initial_energy=a[2]))
 print(rows[-1],flush=True)
out=dict(controls={'width':width,'epsilon':eps,'zeta':zeta,'clock_sigma':.08,'speed':1,'duration':T},rows=rows)
p.with_name('results.json').write_text(json.dumps(out,indent=2)+'\n');print('CP09_BOUNDED_CLOCK_EXECUTION_COMPLETE')
