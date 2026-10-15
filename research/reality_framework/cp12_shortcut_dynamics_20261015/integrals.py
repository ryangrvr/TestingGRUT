"""Unitary propagation and several independent analytic residual upper sums.

For each fixed action L, integrate ||L exp(-i A s)v|| using midpoint Taylor
norms through order five and a sixth-derivative remainder. Ordinary floating
point, not outward-rounded arithmetic. A may be shifted by a real scalar.
"""
from math import ceil, factorial
import json,time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply


def schur_norm(L):
    if not L.nnz:return 0.
    L=abs(L)
    return float(np.max(np.asarray(L.sum(axis=0)))*np.max(np.asarray(L.sum(axis=1))))**.5


class Propagator:
    def __init__(self,H,actions,shift,dt=1/32,batch=4):
        self.A=H-shift*eye(H.shape[0],format='csr');self.generator=-1j*self.A
        self.trace=-1j*complex(self.A.diagonal().sum());self.actions=actions
        self.norms={key:schur_norm(L) for key,L in actions.items()}
        self.dt=dt;self.batch=batch

    def run(self,V,times):
        vector=V.ndim==1
        if vector:V=V[:,None]
        W=V.copy();h2=None
        for order in range(1,7):
            W=self.A@W
            if order==2:h2=np.linalg.norm(W,axis=0)
        h6=np.linalg.norm(W,axis=0)
        errors={key:np.zeros(V.shape[1]) for key in self.actions}
        low={key:np.zeros(V.shape[1]) for key in self.actions}
        clock=0.;answers={};last=time.monotonic()
        for end in times:
            assert end>=clock
            if end>clock:
                steps=max(1,ceil((end-clock)/self.dt));h=(end-clock)/steps
                coefficients=[h**(k+1)/(2**k*factorial(k+1)) for k in range(6)]
                for j in range(0,steps,self.batch):
                    count=min(self.batch,steps-j)
                    grid=expm_multiply(self.generator,V,start=0.,stop=count*h,num=2*count+1,traceA=self.trace)
                    for k in range(count):
                        W=grid[2*k+1]
                        for order,c in enumerate(coefficients):
                            if order:W=self.A@W
                            for key,L in self.actions.items():
                                norm=np.linalg.norm(L@W,axis=0)
                                errors[key]+=c*norm
                                if order<2:low[key]+=c*norm
                        for key,b in self.norms.items():
                            errors[key]+=h**7/(2**6*factorial(7))*b*h6
                            low[key]+=h**3/24*b*h2
                    V=grid[-1].copy()
                    if time.monotonic()-last>15:
                        print(json.dumps(dict(stage='propagation',dimension=self.A.shape[0],
                            time=clock+(j+count)*h,
                            residuals={key:np.minimum(errors[key],low[key]).tolist() for key in self.actions})),flush=True)
                        last=time.monotonic()
            clock=end
            answers[end]=(V[:,0].copy() if vector else V.copy(),
                {key:np.minimum(errors[key],low[key]).copy() for key in self.actions})
        return answers
