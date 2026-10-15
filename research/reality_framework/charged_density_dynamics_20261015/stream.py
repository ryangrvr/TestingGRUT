"""Memory-bounded exponential propagation with full outgoing-error sums."""
from __future__ import annotations
from math import factorial,ceil
import time,json
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply


class StreamEvolution:
    def __init__(self,H,B,shift,max_dt=1/32,batch=4):
        self.A=H-shift*eye(H.shape[0],format='csr')
        self.B=B;self.max_dt=max_dt;self.batch=batch
        self.trace=complex(self.A.diagonal().sum())
        self.generator=-1j*self.A
        if B.shape[0] and B.nnz:
            absolute=abs(B)
            self.bn=float(np.max(np.asarray(absolute.sum(axis=0)))*
                          np.max(np.asarray(absolute.sum(axis=1))))**.5
        else:self.bn=0.

    def run(self,V,times):
        vector=V.ndim==1
        if vector:V=V[:,None]
        start=V.copy();W=V.copy();h2=None
        for k in range(1,7):
            W=self.A@W
            if k==2:h2=np.linalg.norm(W,axis=0)
        h6=np.linalg.norm(W,axis=0)
        error=np.zeros(V.shape[1]);old=np.zeros(V.shape[1]);clock=0.;answers={}
        last_report=time.monotonic()
        for end in times:
            if end==clock:
                answers[end]=(V[:,0].copy() if vector else V.copy(),np.minimum(error,old).copy());continue
            steps=max(1,ceil((end-clock)/self.max_dt));h=(end-clock)/steps
            coefficients=[h**(k+1)/(2**k*factorial(k+1)) for k in range(6)]
            tail=h**7/(2**6*factorial(7))*self.bn*h6
            for j in range(0,steps,self.batch):
                count=min(self.batch,steps-j)
                grid=expm_multiply(self.generator,V,start=0.,stop=count*h,
                    num=2*count+1,traceA=-1j*self.trace)
                for k in range(count):
                    W=grid[2*k+1];local=np.zeros(V.shape[1])
                    for order,c in enumerate(coefficients):
                        if order:W=self.A@W
                        norm=np.linalg.norm(self.B@W,axis=0)
                        local+=c*norm
                        if order<2:old+=c*norm
                    error+=local+tail
                    old+=h**3/24*self.bn*h2
                V=grid[-1].copy()
                if time.monotonic()-last_report>15:
                    print(json.dumps(dict(stage='evolution',dimension=self.A.shape[0],time=clock+(j+count)*h,error=np.minimum(error,old).tolist())),flush=True)
                    last_report=time.monotonic()
            clock=end
            answers[end]=(V[:,0].copy() if vector else V.copy(),np.minimum(error,old).copy())
        return answers


def probability_interval(p,error,allowance=5e-10):
    radius=float(error)+allowance
    return [max(0.,max(0.,float(p))**.5-radius)**2,
            min(1.,(max(0.,float(p))**.5+radius)**2)]
