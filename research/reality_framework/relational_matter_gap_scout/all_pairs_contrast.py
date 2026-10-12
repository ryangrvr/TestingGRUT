#!/usr/bin/env python3
"""Adversarial all-pairs logarithmic phase gas on a FIXED supplied circle.

Comparison isolates propagation locality, not thermal/pressure selection.
The nearest-gap SCOUT law is in check.py. Requires numpy.
"""
import json
from pathlib import Path
from math import pi
import numpy as np

def hessian_log(n,g=1.):
    H=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            weight=g/(4*n*np.sin(pi*(j-i)/n)**2)
            H[i,i]+=weight;H[j,j]+=weight
            H[i,j]-=weight;H[j,i]-=weight
    return H

def hessian_local(n,lam=1.):
    a=2*pi/n
    weight=2*lam/(n*a**3)
    H=np.zeros((n,n))
    for i in range(n):
        j=(i+1)%n
        H[i,i]+=weight;H[j,j]+=weight
        H[i,j]-=weight;H[j,i]-=weight
    return H

def main():
    rows=[]
    for n in (8,16,32,64):
        logeig=np.sort(np.linalg.eigvalsh(hessian_log(n)))
        localeig=np.sort(np.linalg.eigvalsh(hessian_local(n)))
        k=np.arange(n)
        exactlog=np.sort(k*(n-k)/(2*n))
        exactlocal=np.sort(n*n/(pi**3)*np.sin(pi*k/n)**2)
        assert np.max(abs(logeig-exactlog))<1e-10
        assert np.max(abs(localeig-exactlocal))<1e-9
        Lg,Ll=hessian_log(n),hessian_local(n)
        far=n//4
        assert abs(Lg[0,far])>1e-6
        assert abs(Ll[0,far])<1e-12
        rows.append(dict(N=n,
            nearest_phase_spring_omega2_over_omega1=float(np.sqrt(localeig[2]/localeig[1])),
            all_pairs_log_omega2_over_omega1=float(np.sqrt(logeig[2]/logeig[1])),
            distant_pair_hessian=float(Lg[0,far]),
            distant_nearest_gap_hessian=float(Ll[0,far])))
    # In the long-wave region, log all-pairs omega_k~sqrt(k),
    # while gap-local omega_k~k, at fixed circle radius.
    assert rows[-1]['all_pairs_log_omega2_over_omega1']<1.5
    assert rows[-1]['nearest_phase_spring_omega2_over_omega1']>1.99
    print(json.dumps(rows,indent=2))
    print("SCOUT_RG01_ALL_PAIRS_LOCALITY_ADVERSARY_PASS")
if __name__=="__main__":main()
