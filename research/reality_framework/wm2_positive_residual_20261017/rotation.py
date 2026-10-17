"""Finite Fock rotation coefficients enclosed independently of their evaluation."""
from fractions import Fraction as Q
from math import comb,factorial
import numpy as np
import enclosures as E
from layout import S5


def cp(z):return E.point(float(z.real)),E.point(float(z.imag))
def cadd(a,b):return E.iadd(a[0],b[0]),E.iadd(a[1],b[1])
def cneg(a):return E.ineg(a[0]),E.ineg(a[1])
def cmul(a,b):return E.isub(E.imul(a[0],b[0]),E.imul(a[1],b[1])),E.iadd(E.imul(a[0],b[1]),E.imul(a[1],b[0]))
def cscale(a,b):return E.imul(a[0],b),E.imul(a[1],b)
def power(a,n):
    out=cp(1+0j)
    for _ in range(n):out=cmul(out,a)
    return out


def matrices(t,maximum):
    result={};delta=rows=0.0
    for n in range(1,maximum+1):
        mat=S5.fock_two_mode(t,n)
        bound=np.zeros((n+1,n+1))
        for old in range(n+1):
            b=n-old
            for new in range(n+1):
                value=cp(0j)
                for x in range(max(0,new-b),min(old,new)+1):
                    y=new-x
                    term=cmul(cmul(power(cp(t[0,0]),x),power(cp(t[1,0]),old-x)),
                              cmul(power(cp(t[0,1]),y),power(cp(t[1,1]),b-y)))
                    value=cadd(value,cscale(term,E.point(float(comb(old,x)*comb(b,y)))))
                ratio=E.rational_bounds(Q(factorial(new)*factorial(n-new),factorial(old)*factorial(b)))
                scale=(E.sqrt_down(ratio[0]),E.sqrt_up(ratio[1]))
                value=cscale(value,scale)
                z=mat[new,old]
                er=max(E.upper(z.real-value[0][0]),E.upper(value[0][1]-z.real))
                ei=max(E.upper(z.imag-value[1][0]),E.upper(value[1][1]-z.imag))
                bound[new,old]=E.add_up(er,ei)
        delta=max(delta,E.upper(max(np.sum(bound,axis=0).max(),np.sum(bound,axis=1).max())*(1+E.gamma(2*(n+1)))))
        absolute=abs(mat.real)+abs(mat.imag)
        rows=max(rows,E.upper(max(absolute.sum(axis=0).max(),absolute.sum(axis=1).max())*(1+E.gamma(2*(n+1)))))
        result[n]=np.ascontiguousarray(mat)
    # Two-term complex dot, diagonal subtraction and row summation enclosures.
    gram=t.conj().T@t-np.eye(2)
    eta=E.add_up(float(np.max(np.sum(abs(gram.real)+abs(gram.imag),axis=1))),E.gamma(128))
    return result,delta,rows,eta


def fast_matrices(t,maximum):
    """Conservative absolute polynomial evaluation bound, without a density test.

    Each entry contains <=n+1 terms of <=n complex products/powers and exact
    integer binomials. gamma(4096) dominates these operations at n<=5,
    including rounded square roots. The absolute coefficient sum is bounded
    before cancellation; row/column bounds then give an operator norm bound.
    The slower interval implementation independently checks this in controls.
    """
    assert maximum<=5
    result={};delta=rows=0.0
    tmax=E.upper(max(1.0,float(np.max(abs(t.real)+abs(t.imag)))))
    for n in range(1,maximum+1):
        mat=S5.fock_two_mode(t,n)
        absolute=0.0
        for old in range(n+1):
            for new in range(n+1):
                b=n-old
                amount=sum(comb(old,x)*comb(b,new-x) for x in range(max(0,new-b),min(old,new)+1))
                ratio=E.rational_bounds(Q(factorial(new)*factorial(n-new),factorial(old)*factorial(b)))[1]
                absolute=max(absolute,E.upper(amount*E.sqrt_up(ratio)*tmax**n))
        delta=max(delta,E.upper((n+1)*absolute*E.gamma(4096)))
        bound=abs(mat.real)+abs(mat.imag)
        rows=max(rows,E.upper(max(bound.sum(axis=0).max(),bound.sum(axis=1).max())*(1+E.gamma(2*(n+1)))))
        result[n]=np.ascontiguousarray(mat)
    gram=t.conj().T@t-np.eye(2)
    eta=E.add_up(float(np.max(np.sum(abs(gram.real)+abs(gram.imag),axis=1))),E.gamma(128))
    return result,delta,rows,eta
