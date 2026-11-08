"""Frozen finite-volume dictionary and ordered Gaussian reference."""
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'grut_detector_matching_20261107'
sys.path.insert(0,str(OLD))
from continuum import interval_spectrum, all_histories, gaussian_histories, key


def weights(card,n):
    L,s=card['length'],card['sigma']; a=L/(n+1)
    edges=np.arange(1,n+1)*a
    normalization=ndtr(L/(2*s))-ndtr(-L/(2*s))
    v=np.diff(ndtr((edges-L/2)/s))/normalization
    b=1-float(v.sum())
    stagger=(-1.)**np.arange(n-1)
    shift=(b-float(v@stagger))/4
    return a,edges,v,b,shift


def cell_spectrum(card,n,modes):
    L=card['length']; a,edges,v,b,_=weights(card,n)
    k=np.arange(1,modes+1)*np.pi/L
    omega=np.sqrt(k*k+1/np.pi)
    heights=v/a
    fn=np.sqrt(2/L)*(heights@(np.cos(edges[:-1,None]*k)-np.cos(edges[1:,None]*k)))/k
    h2=fn*fn/(2*np.pi*omega)
    mu=1/np.sqrt(np.pi)
    mean=b/4+float(heights@(np.sinh(mu*(edges[1:]-L/2))-np.sinh(mu*(edges[:-1]-L/2))))/(4*mu*np.cosh(mu*L/2))
    B=abs(heights[0])+float(np.abs(np.diff(heights)).sum())+abs(heights[-1])
    return omega,h2,mean,{'boundary_weight':b,'stagger_and_boundary_shift':weights(card,n)[4],
        'weights':v.tolist(),'total_variation_smearing':B,'variance':float(h2.sum()),'mean':mean,
        'tail_variance_bound':B*B*L*L/(2*np.pi**4*modes**2)}


def covariance(omega,h2):
    return {str(t):complex(np.sum(h2*np.exp(-1j*omega*t))) for t in (0,1,2)}


def gaussian_from_cov(card,mean,cov):
    """Auxiliary summaries only; no inferred parameter alters the dictionary."""
    output={}
    for times in card['protocols']:
        n=len(times); c=card['coupling']; words=[]
        def C(dt):
            z=cov[str(int(abs(dt)))];return z if dt>=0 else z.conjugate()
        for bra in itertools.product((0,1),repeat=n):
            for ket in itertools.product((0,1),repeat=n):
                ops=[(i,-c) for i in range(n) if bra[i]]+[(i,c) for i in reversed(range(n)) if ket[i]]
                z=1j*mean*sum(lam for _,lam in ops)
                z+=1j*card['phase']*sum(ket[i]-bra[i] for i in range(n))
                for ai,(i,lam) in enumerate(ops):
                    z-=lam*lam*C(0)/2
                    for j,other in ops[ai+1:]: z-=lam*other*C(times[i]-times[j])
                words.append((bra,ket,np.exp(z)))
        p={}
        for signs in itertools.product((1,-1),repeat=n):
            z=sum(value*math.prod(signs[i]**(bra[i]+ket[i]) for i in range(n)) for bra,ket,value in words)/(4**n)
            if abs(z.imag)>1e-10: raise ValueError('Non-real Gaussian probability')
            p[''.join('+' if s==1 else '-' for s in signs)]=float(z.real)
        output[key(times)]=p
    return output


def tv(p,q):
    if set(p)!=set(q): raise ValueError('Different outcome dictionaries')
    return .5*sum(abs(p[s]-q[s]) for s in p)


def distances(p,q):
    if set(p)!=set(q): raise ValueError('Different histories')
    return {k:tv(p[k],q[k]) for k in p}


def disturbance(p):
    return {k:.5*sum(abs(sum(v for s,v in p[k].items() if s[-1]==outcome)-p['0'][outcome])
                       for outcome in ('+','-')) for k in ('0,1','0,2')}


def references(card,n):
    o,h,_,_,meta=interval_spectrum(card,card['modes'],card['quadrature_epsabs'])
    target,_=all_histories(card,o,h,meta['mean'])
    oc,hc,mc,cm=cell_spectrum(card,n,card['modes'])
    cell,_=all_histories(card,oc,hc,mc)
    of,hf,mf,fm=cell_spectrum(card,n,card['refined_modes'])
    fine,_=all_histories(card,of,hf,mf)
    ot,ht,_,_,mt=interval_spectrum(card,card['refined_modes'],card['quadrature_epsabs'])
    target_fine,_=all_histories(card,ot,ht,mt['mean'])
    r=max(map(len,card['protocols']))
    return {'target_histories':target,'cell_histories':cell,'cell_metadata':cm,'target_metadata':meta,
        'cell_mode_doubling_tv':distances(cell,fine),'target_mode_doubling_tv':distances(target,target_fine),
        'maximum_history_tail_bound':2**r*r*r*card['coupling']**2*max(cm['tail_variance_bound'],meta['tail_variance_bound']),
        'target_disturbance':disturbance(target),'cell_disturbance':disturbance(cell)},covariance(oc,hc)


def verify(card):
    coefficient_error=0.; partition=0.; gap_error=0.
    for n in card['N']:
        a,edges,v,b,shift=weights(card,n)
        partition=max(partition,abs(v.sum()+b-1))
        matrix=np.diag(np.full(n-1,1/(2*a)),1)+np.diag(np.full(n-1,1/(2*a)),-1)
        ev=np.linalg.eigvalsh(matrix)
        gap_error=max(gap_error,abs(ev[n//2]-ev[n//2-1]-2*np.sin(np.pi/(2*(n+1)))/a))
        omega,h2,_,_=cell_spectrum(card,n,12)
        for j in range(1,13):
            k=j*np.pi/card['length']
            direct=sum(vl/a*quad(lambda z:np.sqrt(2/card['length'])*np.sin(k*z),lo,hi,epsabs=1e-13)[0]
                       for lo,hi,vl in zip(edges[:-1],edges[1:],v))
            coefficient_error=max(coefficient_error,abs(direct*direct/(2*np.pi*omega[j-1])-h2[j-1]))
    o=np.array([.7,1.3,2.1]);h=np.array([.08,.03,.005]);mean=.17
    direct,_=all_histories(card,o,h,mean)
    auxiliary=gaussian_from_cov(card,mean,covariance(o,h))
    return {'partition_error':partition,'free_hopping_gap_error':gap_error,
            'cell_coefficient_quadrature_error':coefficient_error,
            'gaussian_word_crosscheck_tv':max(distances(direct,auxiliary).values())}


if __name__=='__main__':
    print(json.dumps(verify(json.loads((HERE/'protocol_b1_01.json').read_text())),indent=2))
