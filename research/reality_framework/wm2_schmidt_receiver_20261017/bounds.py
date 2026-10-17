"""Small positive enclosures under the S8 IEEE/reduction arithmetic model."""
from functools import lru_cache
from fractions import Fraction as Q
from math import nextafter,sqrt,inf
from pathlib import Path
import sys
import numpy as np

S7=Path(__file__).resolve().parent.parent/'wm2_positive_residual_20261017'
sys.path.append(str(S7))
import native
import enclosures as OLD
from layout import S5


def up(x):return nextafter(float(x),inf)
def down(x):return nextafter(float(x),-inf)


@lru_cache(None)
def gamma(n,bits=53):
    t=Q(int(n),2**bits)
    assert t<Q(1,2)
    return OLD.rational_bounds(t/(1-t))[1]


def positive(raw,n,bits=native.LONG_DIGITS):
    g=up(gamma(4*max(1,n)+8,bits)+(gamma(2) if bits>53 else 0))
    pad=(4*max(1,n)+8)*float.fromhex('0x1p-1022')
    return max(0.,down((raw-pad)/(1+g))),up((raw+pad)/(1-g))


def norm2(x):return positive(native.raw_norm(np.asarray(x)),x.size)
def norm(x):return up(sqrt(norm2(x)[1]))


def scalar_norm(values):
    values=np.asarray(values,dtype=float)
    return norm(values.astype(complex))


def gram(x):
    m,r=x.shape
    if not r:return 0.
    raw=x.conj().T@x;raw.flat[::r+1]-=1
    row=float((abs(raw.real)+abs(raw.imag)).sum(axis=1).max())
    col=float((abs(x.real)**2+abs(x.imag)**2).sum(axis=0).max())
    col=up(col*(1+gamma(16*m+16)))
    return up(row*(1+gamma(4*r+8))+r*gamma(8*m+16)*col+gamma(16)*(1+col))


def hermitian_norm(x,construction=0.):
    """Any eigenvector output is checked, not assumed to be exact."""
    x=(x+x.conj().T)/2
    values,vectors=np.linalg.eigh(x)
    eta=gram(vectors)
    if eta>=1e-6:raise RuntimeError('eigenvector Gram certificate failed')
    residual=x@vectors-vectors*values
    n=len(values);a=float(np.max(abs(values)))
    row=float(np.max(np.sum(abs(x.real)+abs(x.imag),axis=1)))
    rounding=gamma(64*n+64)*sqrt(max(1,n))*(1+row+a)
    error=up(norm(residual)+rounding)
    upper=up((a*up(sqrt(1+eta))+error)/down(sqrt(1-eta))+construction)
    return upper,dict(observed_maximum_eigenvalue=a,Gram_upper=eta,
        eigenvector_residual_upper=error,construction_upper=construction,norm_upper=upper)


def physical_probability(interval,error):return OLD.physical_probability(interval,error)
