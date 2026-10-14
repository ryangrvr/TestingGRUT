"""Observable error from forward and backward full-law leakage.

All insertions and observables commute with the electric projection. The
observable need not commute with H, the source, or any density matrix.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_unmasked_formation'))
from run import error_evolve as parent_error_evolve


def error_evolve(H,B,V,end,steps,shift):
    """Same analytic upper sum; reuse one exponential time-grid operation.

    Copies the endpoint so the returned vector does not retain the time grid.
    Grid memory is explicit: (2*steps+1)*dimension*columns complex scalars.
    """
    A=H-shift*eye(H.shape[0],format='csr')
    vector=V.ndim==1
    if vector:V=V[:,None]
    dt=end/steps
    if B.shape[0]:
        bn=np.sqrt(float(np.max(np.asarray(abs(B).sum(axis=0))))*float(np.max(np.asarray(abs(B).sum(axis=1)))))
    else:bn=0.
    h2=np.linalg.norm(A@(A@V),axis=0)
    trace=complex(A.diagonal().sum())
    grid=expm_multiply(-1j*A,V,start=0.,stop=end,num=2*steps+1,traceA=-1j*trace)
    bound=np.zeros(V.shape[1])
    for j in range(steps):
        mid=grid[2*j+1]
        bound+=dt*np.linalg.norm(B@mid,axis=0)+dt*dt/4*np.linalg.norm(B@(A@mid),axis=0)+dt**3/24*bn*h2
    out=grid[-1].copy()
    return (out[:,0],float(bound[0])) if vector else (out,bound)


def response(H,B,psi,source,detectors,wait,lag,steps,shift):
    v,ew=error_evolve(H,B,psi,wait,steps,shift)
    Z,ez=error_evolve(H,B,np.column_stack([v,source*v]),lag,steps,shift)
    x,y=Z.T
    ex=ew+float(ez[0]);ey=ew+float(ez[1])
    global_bound=2*(ex+ey+ex*ey)
    # Batch independent backward probes; each column keeps its own bound.
    DX=np.column_stack([O*x for O in detectors])
    DY=np.column_stack([O*y for O in detectors])
    a,e_lag=error_evolve(-H,B,DX,lag,steps,-shift)
    _,e_wait=error_evolve(-H,B,source[:,None]*a,wait,steps,-shift)
    _,e_dual_x=error_evolve(-H,B,DY,wait+lag,2*steps,-shift)
    rows=[]
    for j,O in enumerate(detectors):
        # Backward path of the source-inserted forward error y.
        dual_y=float(e_lag[j]+e_wait[j])
        # Backward path of the uninserted forward error x.
        dual_x=float(e_dual_x[j])
        # Global detector norm on the UNREGULATED space, not merely its
        # maximum over the retained states. These detector averages are <=1.
        onorm=1.
        bilinear=ey*dual_y+ex*dual_x+onorm*ex*ey
        rows.append(dict(projected_response=float(2*np.vdot(x,O*y).imag),
                         global_absolute_error_upper=global_bound,
                         observable_absolute_error_upper=2*bilinear,
                         forward_x_error_upper=ex,forward_y_error_upper=ey,
                         dual_y_error_upper=dual_y,dual_x_error_upper=dual_x,
                         detector_norm=onorm))
    return rows,Z,v
