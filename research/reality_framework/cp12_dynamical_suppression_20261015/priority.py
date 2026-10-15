"""Deterministic outgoing-fluence proxy used ONLY for numerical basis selection."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
from boundary import RF,load_model
from operators import flow


def main():
    a=argparse.ArgumentParser();a.add_argument('--n',type=int,default=6);a.add_argument('--q',type=int,default=12)
    a.add_argument('--budget',type=int,default=12);a.add_argument('--model');a.add_argument('--wait',type=float,default=1.)
    a.add_argument('--lag',type=float,default=.5);a.add_argument('--out',required=True);args=a.parse_args()
    m=load_model(args.n,args.q,args.budget,args.model);H=m['H'];B=m['B'];A=-1j*(H-args.q/2*eye(H.shape[0],format='csr'))
    grid=expm_multiply(A,m['psi'],start=0,stop=args.wait,num=5,traceA=A.diagonal().sum())
    P=np.array([flow(args.n,s)[0].bit_count()>=int(np.ceil(args.q/4)) and
        (s[1]&((1<<args.n)-1)).bit_count()>=int(np.ceil(args.n/4)) for s in m['states']])
    X=P*grid[-1];p=float(np.linalg.norm(X)**2)
    charged=expm_multiply(A,X,start=0,stop=args.lag,num=3,traceA=A.diagonal().sum())
    score=np.sum(abs(B@grid.T)**2,axis=1)+np.sum(abs(B@charged.T)**2,axis=1)/p
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    np.savez(args.out,priority=score,wait=grid[-1],charged=X)
    cum=np.cumsum(np.sort(score)[::-1]);cum/=cum[-1]
    print(json.dumps(dict(out=args.out,formation_probability=p,
        counts={str(a):int(np.searchsorted(cum,a)+1) for a in [.9,.99,.999,.9999]})))


if __name__=='__main__':main()
