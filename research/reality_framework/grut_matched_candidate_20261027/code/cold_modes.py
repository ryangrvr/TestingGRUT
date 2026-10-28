"""Exact coupled cold-fluid dispersion for the changed potential.

The density sweep is an auxiliary mechanism test, not a measured onset population.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial import polynomial as poly

ROOT=Path(__file__).resolve().parents[1]


def main():
    a,b,c=.1,3.,.002/.01
    rows=[]
    for F in (.1,.5,1.,1.5,2.,2.63):
        S=brentq(lambda s:b*(s*s-1)*s-3*c*(s-1)**2+F*F*s,0.,1.)
        omega=np.sqrt(S*S+a*F*F)
        D=b*(3*S*S-1)-6*c*(S-1)+F*F
        ks=np.linspace(0,3,1201);growth=[];residual=0.
        for k in ks:
            k2=k*k
            # [(z+k²+2aF²)(z+k²)+4Omega² z](z+k²+D)-4F²S²(z+k²)=0.
            p=poly.polymul([k2+2*a*F*F,1.],[k2,1.]);p[1]+=4*omega*omega
            coeff=poly.polymul(p,[k2+D,1.]);coeff[:2]-=4*F*F*S*S*np.array([k2,1.])
            z=poly.polyroots(coeff)
            frequencies=np.sqrt(z.astype(complex))
            dominant=frequencies[np.argmax(frequencies.real)]
            value=max(0.,float(dominant.real))
            growth.append(value)
            if value>0:
                nu=dominant
                matrix=np.array([[nu*nu+k2+2*a*F*F,-2*omega*nu,2*F*S],
                                 [2*omega*nu,nu*nu+k2,0.],
                                 [2*F*S,0.,nu*nu+k2+D]])
                residual=max(residual,abs(np.linalg.det(matrix))/max(1,np.linalg.norm(matrix)**3))
        peak=int(np.argmax(growth));kp=ks[peak]
        rows.append({'supplied_F':F,'self_consistent_S':S,'Omega':omega,'neutral_curvature_D':D,
                     'cold_longwave_instability_criterion':bool(a*D<2*S*S),
                     'peak_growth_over_mX':float(growth[peak]),'peak_k_over_mX':float(kp),
                     'wavelength_cube_charge_diagnostic':float(omega*F*F/.01*(2*np.pi/kp)**3) if kp else None,
                     'determinant_residual':float(residual)})
    (ROOT/'COLD_MODE_RESULTS.json').write_text(json.dumps({'scope':__doc__,'rows':rows,'formed_population':None},indent=2)+'\n')
    print(json.dumps({'cold_mode_rows':rows}))


if __name__=='__main__':main()
