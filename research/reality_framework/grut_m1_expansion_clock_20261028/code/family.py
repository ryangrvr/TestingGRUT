"""Local finite-Q stationary energy spline and held-out first-law checks."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from mass_analysis import stationary

ROOT=Path(__file__).resolve().parents[1]
def main():
    q=np.array([180000.,185000.,190000.,195000.,200000.,205000.])
    rows=[stationary(float(x)) for x in q]
    es=np.array([r['energy_mX'] for r in rows])
    os=np.array([r['Omega'] for r in rows])
    spline=CubicHermiteSpline(q,es,os)
    holds=[]
    for h in (182500.,187111.,192500.,197500.,201600.,202500.):
        r=stationary(h)
        r['energy_interpolation_error_mX']=float(spline(h)-r['energy_mX'])
        r['omega_interpolation_error']=float(spline(h,1)-r['Omega'])
        holds.append(r)
    out={'nodes':rows,'holdouts':holds,'method':'Cubic Hermite interpolation of stationary energy using independently solved frequencies as derivatives',
         'scope':'stationary energy and frequency only; scattering strengths remain tied to the fixed reference profile'}
    (ROOT/'FAMILY.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
