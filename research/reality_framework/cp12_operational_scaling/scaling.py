"""Full frozen CP12 versus compressed graph block and bare conditional dynamics."""
from pathlib import Path
import sys,itertools,json,time,argparse
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,hop,GAMMA

def run(n,L,long=False):
 start=time.time();q=n;edges=list(itertools.combinations(range(n),2));M=len(edges);eid={e:i for i,e in enumerate(edges)}
 ring=sum(1<<eid[tuple(sorted((i,(i+1)%n)))] for i in range(n))
 hostile=sum(1<<eid[0,i] for i in range(1,n)) | (1<<eid[1,2])
 H,states,lookup,degree=build(n,q,L);masks=np.array([s[0] for s in states]);bits=np.array([s[1] for s in states]);boundary=np.array([max(abs(e) for e in s[2])==L for s in states]);B=((bits>>0)&1).astype(float)
 detectors=[1,n//2];A=np.array([((bits>>j)&1).astype(float) for j in detectors]);b2=GAMMA**2*q*(M-q)/M**2;rows=[]
 for kind,mask in [('ring',ring),('hub',hostile)]:
  ids=np.flatnonzero(masks==mask);block=H[ids][:,ids];bare=block-diags(GAMMA/M*degree[ids]);sea=sum(1<<(n+e) for e in range(M));ef=(0,)*(2*M)
  vacuum=lookup[mask,sea,ef];nb,sg=hop(sea,0,n+eid[0,1]);fields=list(ef);fields[2*eid[0,1]]=1;pair=lookup[mask,nb,tuple(fields)]
  psi=np.zeros(len(states),complex);psi[vacuum]=1/np.sqrt(2);psi[pair]=1j*sg/np.sqrt(2)
  startvec=np.column_stack((psi,B*psi));full=startvec.copy();fixed=startvec[ids].copy();plain=fixed.copy();last=0;energy=np.vdot(psi,H@psi).real
  for T in ([.1,.4,1.,float(np.sqrt(n)),float(2*np.sqrt(n)),float(np.sqrt(n)/GAMMA)] if long else [.1,.4,1.,float(np.sqrt(n)),float(2*np.sqrt(n))]):
   dt=T-last;full=expm_multiply(-1j*dt*H,full);fixed=expm_multiply(-1j*dt*block,fixed);plain=expm_multiply(-1j*dt*bare,plain);last=T
   assert abs(np.linalg.norm(full[:,0])-1)<1e-9
   p=abs(full[:,0])**2;leak=float(1-p[ids].sum());assert leak<=b2*T*T+1e-9
   assert abs(np.vdot(full[:,0],H@full[:,0]).real-energy)<1e-8
   record=dict(N=n,L=L,preparation=kind,T=T,leakage=leak,leakage_bound=b2*T*T,boundary=float(p@boundary),detectors=[])
   for j,a in zip(detectors,A):
    response=lambda v,d:float(-2*np.imag(np.vdot(v[:,0],d*v[:,1])))
    rf=response(full,a);rb=response(fixed,a[ids]);rp=response(plain,a[ids]);efull=float(p@a);eblock=float(abs(fixed[:,0])**2@a[ids])
    assert abs(rf-rb)<=4*b2*T*T+1e-8
    assert abs(efull-eblock)<=2*b2*T*T+1e-8
    record['detectors'].append(dict(vertex=j,response_full=rf,response_block=rb,response_bare=rp,response_error_bound=4*b2*T*T,expectation_full=efull,expectation_block=eblock))
   rows.append(record);print(n,L,kind,T,leak,record['boundary'],flush=True)
 result=dict(N=n,Q=q,L=L,dimension=len(states),seconds=time.time()-start,rows=rows)
 Path(__file__).with_name(f'scaling_N{n}_L{L}'+('_long' if long else '')+'.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4);p.add_argument('--L',type=int,default=2);p.add_argument('--long',action='store_true');a=p.parse_args();run(a.n,a.L,a.long)
