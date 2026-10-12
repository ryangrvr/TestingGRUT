#!/usr/bin/env python3
"""Conditional finite-cycle spectrum with massive incidence fermions.

A deliberately counterfactual mass check; frozen CP-11 had m=1, new
SCOUT-RG02 theorem uses m=0. Does NOT assert a complete positive-mass
graph-quantum ground theorem or include finite-temperature excited fermions.
"""
import math,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad

def excess(l,m,t=1):
    # bulk trace of sqrt(m^2+L) minus discrete l-cycle trace
    f=lambda q: math.sqrt(m*m+4*t*t*math.sin(q/2)**2)
    bulk=quad(f,0,2*math.pi,epsabs=1e-12,epsrel=1e-12)[0]/(2*math.pi)
    q=2*math.pi*np.arange(l)/l
    return l*bulk-sum(f(float(x)) for x in q)
def main():
    rows=[]
    for m in (0.,.25,1.,2.):
        vals=[(l,excess(l,m)) for l in (1,2,3,4,5,8,12)]
        assert all(v>0 for l,v in vals)
        assert all(vals[i][1]>vals[i+1][1] for i in range(len(vals)-1))
        rows.append(dict(mass=m,excess={str(l):v for l,v in vals}))
    assert abs(rows[2]['excess']['8']-.000143949529589)<1e-10
    result=dict(controls=rows,
        general_analytic_claim=("For positive mass, square-root-Laplacian Fourier "
          "coefficients away from zero are strictly negative; positive "
          "cycle excess is exponentially small as cycle length grows."),
        warning="Numerical monotonicity is checked, not a universal m>0 monotonicity theorem.")
    Path(__file__).with_name("mass_results.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print("SCOUT_RG02_MASS_ROBUSTNESS_PASS")
if __name__=="__main__":main()
