#!/usr/bin/env python3
"""Weighted Birman-Schwinger angular-tail bound; see SPECTRAL_NOTES.md."""
import hashlib
import json
from pathlib import Path
import numpy as np
import spectral as sp
from audit import interval, product, square

HERE=Path(__file__).resolve().parent


def weighted_enclosure(R):
    edges=np.r_[sp.SPL.x[(sp.SPL.x>=0)&(sp.SPL.x<R)],R]
    bound=0.; at=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        h=hi-lo
        d=sp.SPL(lo)[[0,2]]; c=sp.SPL(lo,1)[[0,2]]
        b=sp.SPL(lo,2)[[0,2]]/2; a=sp.SPL(lo,3)[[0,2]]/6
        bern=np.array([d,d+c*h/3,d+2*c*h/3+b*h*h/3,d+c*h+b*h*h+a*h**3])
        F,S=[(float(bern[:,j].min())-1e-10,float(bern[:,j].max())+1e-10) for j in range(2)]
        f2,s2=square(F),square(S)
        A=interval(s2,3,f2,offset=-3); B=interval(s2,1,f2,.2,-1)
        C=(.1*f2[0],.1*f2[1]); D=tuple(2*x for x in product(F,S)); E=(.01*(s2[0]-1),.01*(s2[1]-1))
        diag=[A,interval(B,1,C),interval(B,1,C,-1),E]
        mid=np.diag([(v[0]+v[1])/2 for v in diag]); mid[0,1]=mid[1,0]=(D[0]+D[1])/2
        rad=np.diag([(v[1]-v[0])/2 for v in diag]); rad[0,1]=rad[1,0]=(D[1]-D[0])/2
        val=hi**2*(max(abs(np.linalg.eigvalsh(mid)))+max(np.sum(rad,axis=1)))
        if val>bound:bound=float(val);at=float(hi)
    return {'R':R,'Wmax_upper_spline':bound,'cell_upper_radius_at_max':at}


def main():
    raw=(HERE/'tail_refinement_protocol.json').read_bytes(); protocol=json.loads(raw)
    enclosures=[weighted_enclosure(R) for R in [20.,25.,50.]]
    W=protocol['Wmax_used']; J=protocol['Irad3_used']; Z=protocol['frequency_cutoff']; L=protocol['angular_cutoff']
    ell=np.arange(L+1,100001,dtype=float); ratio=W/(ell*(ell+1))
    value=Z*J/(3*np.pi)*sum(ratio/((1-ratio)*(2*ell+3)*(2*ell-1)))
    M=100000.; rr=W/M**2
    tail=Z*J/(3*np.pi)*W/(9*(1-rr)*M**3)
    result={'protocol_sha256':hashlib.sha256(raw).hexdigest(),'enclosures':enclosures,'angular_bound':float(value+tail),'series_tail_bound':float(tail),'first_omitted_norm_ratio':float(ratio[0]),'checks':{'Wmax_enclosed':all(x['Wmax_upper_spline']<W for x in enclosures),'norm_ratio_below_one':max(ratio)<1},'scope':'Analytical bound for supplied spline, conditional on padded floating-point enclosures. Not an exact continuum-profile certificate.'}
    result['checks']={k:bool(v) for k,v in result['checks'].items()};result['all_checks_pass']=all(result['checks'].values())
    (HERE/'tail_refinement_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
