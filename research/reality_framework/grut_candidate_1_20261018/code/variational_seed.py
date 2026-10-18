#!/usr/bin/env python3
"""Exploratory initial-profile calculation, not an Euler-Lagrange solution.

Fixed a=1,b=.1,Qbar=500. This calculation was used to suggest an initial
guess to the boundary-value solver. It is not a physical coupling search.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize

r=np.linspace(0,40,4001); weight=4*np.pi*r*r; Q=500.
def energy(x,full=False):
    A,R,B=x
    z=np.exp(-.5*(r/R)**2);F=A*z;S=1-B*z
    Fp=-r/R**2*F;Sp=B*r/R**2*z
    I=simpson(weight*F**2,x=r)
    U=.25*(S*S-1)**2+.5*S*S*F*F+.025*F**4
    E=Q*Q/(2*I)+simpson(weight*(.5*(Fp*Fp+Sp*Sp)+U),x=r)
    return (float(E),float(Q/I)) if full else E

if __name__=='__main__':
    opt=minimize(energy,[2.,3.,.99],bounds=[(.1,3.2),(.3,12),(.01,1)],
                 method='Nelder-Mead',options={'maxiter':2000,'xatol':1e-9,'fatol':1e-9})
    E,omega=energy(opt.x,True)
    result={'status':'exploratory_trial_energy_only','a':1.,'b':.1,'Qbar':Q,
            'parameters_A_R_B':opt.x.tolist(),'Ebar':E,'Omega':omega,
            'Ebar_over_Qbar':E/Q,'optimizer_success':bool(opt.success),
            'qualification':'A finite-energy trial below the free-charge threshold; no stationarity or stability certificate.'}
    Path(__file__).with_name('variational_seed_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
