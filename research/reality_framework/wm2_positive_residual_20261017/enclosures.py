"""Outward scalar/positive bounds under the declared IEEE arithmetic model."""
from fractions import Fraction as Q
from math import factorial,sqrt
import numpy as np
import native

U=Q(1,2**53)
UF=float.fromhex('0x1p-1022')
ETA=1e-9


def up(x):return np.nextafter(x,np.inf)
def down(x):return np.nextafter(x,-np.inf)
def upper(x):return float(up(x))
def lower(x):return float(down(x))
def add_up(a,b):return upper(float(a)+float(b))
def mul_up(a,b):return upper(float(a)*float(b))


def rational_bounds(q):
    q=Q(q);a=b=float(q)
    while Q.from_float(a)>q:a=float(down(a))
    while Q.from_float(b)<q:b=float(up(b))
    return a,b


def sqrt_bounds(x):
    x=max(0.0,float(x));q=Q.from_float(x)
    a=b=sqrt(x)
    while Q.from_float(a)**2>q:a=float(down(a))
    while Q.from_float(b)**2<q:b=float(up(b))
    return a,b


def sqrt_up(x):return sqrt_bounds(x)[1]
def sqrt_down(x):return sqrt_bounds(x)[0]
def gamma(n,bits=53):
    t=Q(int(n),2**bits)
    assert t<Q(1,2)
    return rational_bounds(t/(1-t))[1]


def positive_raw_bounds(raw,n,bits=53,operations=4):
    g=gamma(operations*max(1,int(n))+8,bits)
    # Account separately for conversion of extended-precision sums to binary64.
    if bits>53:g=add_up(g,gamma(2))
    pad=mul_up(operations*max(1,n)+8,UF)
    return max(0.0,lower((raw-pad)/(1+g))),upper((raw+pad)/(1-g))


def norm2_bounds(x):
    return positive_raw_bounds(native.raw_norm(x),x.size,native.LONG_DIGITS)


def norm_up(x):return sqrt_up(norm2_bounds(x)[1])
def norm_interval(x):
    lo,hi=norm2_bounds(x);return sqrt_down(lo),sqrt_up(hi)


def iadd(a,b):return down(a[0]+b[0]),up(a[1]+b[1])
def ineg(a):return -a[1],-a[0]
def isub(a,b):return iadd(a,ineg(b))
def imul(a,b):
    p=np.stack(np.broadcast_arrays(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]))
    return down(p.min(axis=0)),up(p.max(axis=0))
def idiv(a,b):
    assert np.all(np.asarray(b[0])>0)
    return imul(a,(down(1/b[1]),up(1/b[0])))
def point(a):return np.asarray(a),np.asarray(a)


def pi_interval():
    def atan_interval(base,n):
        s=sum((Q((-1)**j,(2*j+1)*base**(2*j+1)) for j in range(n)),Q(0))
        t=Q((-1)**n,(2*n+1)*base**(2*n+1))
        return min(s,s+t),max(s,s+t)
    a,b=atan_interval(5,40),atan_interval(239,12)
    lo,hi=16*a[0]-4*b[1],16*a[1]-4*b[0]
    return rational_bounds(lo)[0],rational_bounds(hi)[1]


def trig_interval(x,cosine=False):
    square=imul(x,x)
    terms=16
    def coefficient(j):return rational_bounds(Q((-1)**j,factorial(2*j+(0 if cosine else 1))))
    p=coefficient(terms-1)
    for j in range(terms-2,-1,-1):p=iadd(imul(p,square),coefficient(j))
    if not cosine:p=imul(p,x)
    xmax=float(np.max(np.maximum(abs(x[0]),abs(x[1]))))
    bound=1.0
    for _ in range(32):bound=mul_up(bound,xmax)
    bound=upper(bound/factorial(32))
    return down(p[0]-bound),up(p[1]+bound)


def center(interval):
    mid=np.asarray(interval[0])+(np.asarray(interval[1])-interval[0])/2
    radius=np.maximum(up(mid-interval[0]),up(interval[1]-mid))
    return mid,radius


def receiver(modes):
    angle=imul(pi_interval(),idiv(point(np.arange(1,modes+1,dtype=float)),point(float(modes+1))))
    freq=imul(point(-2.0),trig_interval(angle,True))
    scale=rational_bounds(Q(2,modes+1))
    scale=(sqrt_down(scale[0]),sqrt_up(scale[1]))
    contact_interval=imul(scale,trig_interval(angle))
    contact,radius=center(contact_interval)
    return contact.astype(float),norm_up(radius.astype(complex)),freq


def phases(energy_interval,t):
    x=imul(energy_interval,point(float(t)))
    co,cr=center(trig_interval(x,True));si,sr=center(trig_interval(x,False))
    return np.asarray(co-1j*si,dtype=complex),float(np.max(up(cr+sr)))


def gram_upper(w):
    r=w.shape[1]
    if r==0:return 0.0
    g=w.conj().T@w-np.eye(r)
    raw=np.sum(abs(g.real)+abs(g.imag),axis=1)
    rows=max(positive_raw_bounds(float(v),2*r)[1] for v in raw)
    col=max(norm2_bounds(w[:,j])[1] for j in range(r))
    dot=mul_up(r,mul_up(gamma(8*w.shape[0]+8),col))
    return add_up(rows,add_up(dot,mul_up(r*w.shape[0]*16,UF)))


def contact_bounds(w,z,c,c_error,eta):
    rank=w.shape[1]
    if not rank:return add_up(norm_up(c.astype(complex)),c_error),c_error
    residual=c-w@z
    zn=norm_up(z);cn=norm_up(c.astype(complex))
    # Positive norm of the directly formed vector plus dot/subtraction bounds.
    dot=mul_up(gamma(8*rank+8),mul_up(max(1,rank)*(1+eta),zn))
    subtraction=mul_up(gamma(4),add_up(cn,mul_up(1+eta,zn)))
    sigma=add_up(norm_up(residual),add_up(c_error,add_up(dot,subtraction)))
    zround=mul_up(gamma(8*w.shape[0]+8),max(1,rank)*(1+eta)*cn)
    projected=add_up(eta,add_up(zround,mul_up(1+eta,c_error)))
    return sigma,projected


def probability_interval(raw_event,raw_norm,n):
    a,b=positive_raw_bounds(raw_event,n,native.LONG_DIGITS)
    c,d=positive_raw_bounds(raw_norm,n,native.LONG_DIGITS)
    return max(0.0,lower(a/d)),min(1.0,upper(b/c))


def physical_probability(interval,e):
    lo,hi=interval
    rev=max(0.0,lower(sqrt_down(lo)-e));rev=max(0.0,lower(rev*rev))
    other=upper(sqrt_up(max(0.0,upper(1-lo)))+e)
    low=max(0.0,rev,lower(1-upper(other*other)))
    a=upper(sqrt_up(hi)+e);upper1=min(1.0,upper(a*a))
    b=max(0.0,lower(sqrt_down(max(0.0,lower(1-hi)))-e))
    upper2=min(1.0,upper(1-lower(b*b)))
    return low,min(upper1,upper2)
