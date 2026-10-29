"""Independent finite-width versus NWA energy-integral checks at fixed angle."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from local_rates import GH,MH,GAMMA

ROOT=Path(__file__).resolve().parents[1]

def calculate(Ed,mu):
    pd=np.sqrt(Ed*Ed-1)
    def sf(E):return 2+2*(Ed*E-pd*np.sqrt(max(E*E-1,0))*mu)
    def fun(E):return sf(E)-6
    # Dense bracketing is independent of the algebraic root formula.
    mesh=np.linspace(1,8,30001);f=np.array([fun(e) for e in mesh]);roots=[]
    for i in range(len(mesh)-1):
        if f[i]*f[i+1]<0:roots.append(brentq(fun,mesh[i],mesh[i+1],xtol=1e-14))
    def source(E):return np.exp(-((E-2.)/1.2)**2)
    def integrand(E):
        p=np.sqrt(max(E*E-1,0));s=sf(E);beta=np.sqrt(max(0,1-4/s));v=p/E
        if v==0:return 0.
        return GH**4*beta/(32*np.pi*E*Ed*((s-6)**2+6*GAMMA**2))*source(E)/v
    knots=[1.,8.];nwa=0
    for E in roots:
        p=np.sqrt(E*E-1);v=p/E;deriv=2*(Ed-pd*mu*E/p)
        nwa+=np.pi*GH**2/(2*Ed*E*abs(deriv))*source(E)/v
        width=MH*GAMMA/abs(deriv)
        knots.extend(x for x in E+width*np.array([-100,-10,-1,0,1,10,100]) if 1<x<8)
    knots=sorted(set(knots))
    exact=sum(quad(integrand,a,b,epsabs=1e-11,epsrel=1e-8,limit=200)[0] for a,b in zip(knots[:-1],knots[1:]))
    return {'Ed':Ed,'relative_cosine':mu,'roots':roots,'finite_width_integral':exact,'nwa_integral':nwa,
            'finite_over_nwa':exact/nwa if nwa else None}

def main():
    rows=[calculate(E,mu) for E,mu in [(1.1,-.8),(1.4,-.5),(1.4,0.),(1.4,.5),(2.5,.9)]]
    out={'rows':rows,'scope':'Fixed-angle smooth synthetic primary spectrum, full s-channel Breit-Wigner vs NWA. This is not a bound on the actual time-dependent daughter-plume calculation or its angular caustics.'}
    (ROOT/'GENERAL_WIDTH_CONTROL.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
