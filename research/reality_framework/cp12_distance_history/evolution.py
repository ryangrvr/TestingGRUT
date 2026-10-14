"""Finite-time exponential propagation and analytic leakage upper sums.

Fifth-order midpoint Taylor bounds the LEAKAGE VECTOR, not the physical
observable's time evolution. Matrix exponentials supply the evolving state.
"""
from math import factorial
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply


class Evolution:
    def __init__(self,H,B,shift,steps=16):
        self.H=H;self.B=B;self.A=H-shift*eye(H.shape[0],format='csr');self.steps=steps
        self.trace=complex(self.A.diagonal().sum())
        self.bn=(float(np.max(np.asarray(abs(B).sum(axis=0))))*float(np.max(np.asarray(abs(B).sum(axis=1)))))**.5 if B.shape[0] else 0.

    def __call__(self,V,end,back=False,groups=None):
        vector=V.ndim==1
        if vector:V=V[:,None]
        if end==0:return (V[:,0].copy() if vector else V.copy()),np.zeros(V.shape[1]),(np.zeros(len(groups)) if groups is not None else 0.)
        selected=[np.asarray(x,dtype=int) for x in groups] if groups is not None else [np.arange(V.shape[1])]
        h=end/self.steps;W=V
        h2=None
        for k in range(1,7):
            W=self.A@W
            if k==2:h2=np.linalg.norm(W,axis=0)
        h6=np.linalg.norm(W,axis=0)
        grid=expm_multiply((1j if back else -1j)*self.A,V,start=0.,stop=end,num=2*self.steps+1,traceA=(1j if back else -1j)*self.trace)
        bound=np.zeros(V.shape[1]);global_bound=np.zeros(len(selected));old=np.zeros(V.shape[1]);old_global=np.zeros(len(selected))
        coefficients=[h**(k+1)/(2**k*factorial(k+1)) for k in range(6)]
        tail=h**7/(2**6*factorial(7))*self.bn
        for j in range(self.steps):
            W=grid[2*j+1];local=np.zeros(V.shape[1]);g=np.zeros(len(selected))
            for k,c in enumerate(coefficients):
                if k:W=self.A@W
                Z=self.B@W;norm=np.linalg.norm(Z,axis=0);total=np.array([np.linalg.norm(norm[t]) for t in selected])
                local+=c*norm;g+=c*total
                if k<2:old+=c*norm;old_global+=c*total
            bound+=local+tail*h6;global_bound+=g+tail*np.array([np.linalg.norm(h6[t]) for t in selected])
            old+=h**3/24*self.bn*h2;old_global+=h**3/24*self.bn*np.array([np.linalg.norm(h2[t]) for t in selected])
        # Both upper sums are valid; take the smaller independently.
        bound=np.minimum(bound,old);global_bound=np.minimum(global_bound,old_global)
        out=grid[-1].copy()
        return (out[:,0] if vector else out),bound,(global_bound if groups is not None else float(global_bound[0]))


def response_bound(ex,ey,dx,dy):
    return 2*(ey*dy+ex*dx+ex*ey)


def dual(E,X,Y,observable,adjoint_x,adjoint_y,wait,lag,ex,ey,pre=None):
    """Any contraction at waiting time, including a recording isometry.

    adjoint_x/y map projected post-wait path vectors to the pre-wait space.
    X,Y may be purification matrices; their inner product is Frobenius.
    """
    vector=X.ndim==1
    DX=observable*X if vector else observable[:,None]*X
    DY=observable*Y if vector else observable[:,None]*Y
    width=1 if vector else X.shape[1]
    Z,_,combined=E(np.column_stack([DX,DY]),lag,True,[range(width),range(width,2*width)])
    zx=Z[:,0] if vector else Z[:,:width]
    zy=Z[:,1] if vector else Z[:,width:]
    evolution_pre=E if pre is None else pre
    _,pre_errors,_=evolution_pre(np.column_stack([adjoint_y(zx),adjoint_x(zy)]),wait,True)
    dyl,dxl=combined;dyw,dxw=pre_errors
    dy=dyl+dyw;dx=dxl+dxw
    value=float(2*np.vdot(X,DY).imag)
    return dict(value=value,error=response_bound(ex,ey,dx,dy),forward_x=ex,forward_y=ey,dual_x=dx,dual_y=dy,interval=[value-response_bound(ex,ey,dx,dy),value+response_bound(ex,ey,dx,dy)])
