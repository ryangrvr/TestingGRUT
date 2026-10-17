"""Direct positive residuals and unknown physical-prefix Gram correction."""
from math import sqrt
import numpy as np
import bounds as E


def add(a,b):return E.up(a+b)
def mul(a,b):return E.up(a*b)
def div(a,b):return E.up(a/E.down(b))


def residual_enclosure(X,K,Z,eta,sum_eta,environment):
    """Exact spans of stored K define the projector; all input arrays are floats.

    BLAS complex reductions obey the declared absolute gamma bounds. Positive
    norms enclose stored arrays; products/subtraction and prefix conditioning
    are additional. This routine assumes neither an exact SVD nor H=I.
    """
    assert 0<=sum_eta<1e-6 and 0<=eta<1e-8
    n,r=K.shape
    if r:
        R=X-K@Z
        kf,zf,xf=E.norm(K),E.norm(Z),E.norm(X)
        scale=add(xf,mul(kf,zf))
        underflow=(32*r+64)*max(1,X.size)*float.fromhex('0x1p-1022')
        rr=add(mul(E.gamma(16*r+32),scale),underflow)
    else:R=X;rr=0.;kf=0.
    rn=E.norm(R);rtrue=add(rn,rr)
    loss=mul(environment,rtrue)
    if not r:return dict(orthogonal_loss_upper=loss,in_range_defect_upper=0.,residual_rounding_upper=rr)
    w0=K.conj().T@R;wn=E.norm(w0)
    kn=E.up(sqrt(add(1.,eta)))
    wr=mul(E.gamma(16*n+32),mul(kf,rn))
    wr=add(wr,(32*n+64)*max(1,w0.size)*float.fromhex('0x1p-1022'))
    wr=add(wr,mul(kn,rr))
    h=div(sum_eta,E.down(1-sum_eta))
    lower=E.down(E.down(1-h)*E.down(1-eta))
    upper=mul(add(1.,h),add(1.,eta))
    condition=div(E.up(sqrt(upper)),lower)
    wtrue=add(add(wn,wr),mul(mul(h,kn),rtrue))
    defect=mul(mul(environment,condition),wtrue)
    return dict(orthogonal_loss_upper=loss,in_range_defect_upper=defect,
        residual_rounding_upper=rr,contraction_rounding_upper=wr,
        unknown_Gram_upper=h,Gram_inverse_factor=condition)
