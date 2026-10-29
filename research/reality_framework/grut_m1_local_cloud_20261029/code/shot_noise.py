"""Poisson-center shot noise of collisionless direct primary clouds.

For homogeneous independent Poisson core centers and identical synchronized
sources, Campbell's theorem gives mean and covariance from a single radial
number-density profile. h daughters, exclusion/formation correlations,
interactions and cosmic expansion are excluded. The pointwise covariance is
not itself a scattering rate, because pair angles and correlations matter.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson
from local_rates import pulse_series

ROOT=Path(__file__).resolve().parents[1]
NCORE=1.2141071422639381e-15


def profile(t_now,r,t,e,F,v,r0):
    """Flux F in counts/time/bin -> local density in counts/m_X^-3."""
    rt=t_now-(r[:,None]-r0)/v[None,:];dt=t[1]-t[0]
    index=np.floor(rt/dt).astype(int);f=rt/dt-index
    valid=(index>=0)&(index<len(t)-1)
    index=np.clip(index,0,len(t)-2);cols=np.arange(len(e))[None,:]
    source=((1-f)*F[index,cols]+f*F[index+1,cols])*valid
    return np.sum(source/v[None,:],axis=1)/(4*np.pi*r*r)


def run(z,times,dr,dt,chunk):
    t,e,plus,minus=pulse_series(z,dt);de=float(np.diff(z['energy_edges'])[0]);plus*=de
    v=np.sqrt(e*e-1)/e;r0=float(z['probe_radius_mX_inverse']);rows=[]
    for T in times:
        # Midpoint radial quadrature resolves narrow travelling source shells.
        nr=int(np.ceil(T/dr));r=r0+(np.arange(nr)+.5)*dr
        one=np.empty(nr);two=np.empty(nr)
        for start in range(0,nr,chunk):
            end=min(start+chunk,nr)
            one[start:end]=profile(T,r[start:end],t,e,plus,v,r0)
            two[start:end]=profile(T,r[start:end],t,e,minus,v,r0)
        vol=4*np.pi*r*r*dr
        I1=float(vol@one);I2=float(vol@two)
        V1=float(vol@(one*one));V2=float(vol@(two*two));cross=float(vol@(one*two))
        target1=float(simpson(plus.sum(axis=1),x=t));target2=float(simpson(minus.sum(axis=1),x=t))
        mean1=NCORE*I1;mean2=NCORE*I2
        rows.append({'mX_t':T,'radial_grid_points':nr,'plus_number_integrated':I1,'minus_number_integrated':I2,
                     'plus_source_number':target1,'minus_source_number':target2,
                     'number_relative_error_plus':I1/target1-1,'number_relative_error_minus':I2/target2-1,
                     'mean_plus_density_over_mX3':mean1,'mean_minus_density_over_mX3':mean2,
                     'poisson_cv_plus':float(np.sqrt(NCORE*V1)/mean1),
                     'poisson_cv_minus':float(np.sqrt(NCORE*V2)/mean2),
                     'same_source_cross_moment_ratio':float(NCORE*cross/(mean1*mean2)),
                     'poisson_mean_product_over_product_of_means':float(1+NCORE*cross/(mean1*mean2))})
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='source_t400');p.add_argument('--dr',type=float,default=20.)
    p.add_argument('--dt',type=float,default=10.);p.add_argument('--chunk',type=int,default=1000);p.add_argument('--tag',required=True);a=p.parse_args()
    z=np.load(ROOT/(a.source+'.npz'));rows=run(z,[10000.,50000.,100000.,200000.,500000.,1000000.],a.dr,a.dt,a.chunk)
    out={'configuration':vars(a),'ncore_over_mX3':NCORE,'rows':rows,
         'scope':'Direct primary radial clouds and independent Poisson core centers only; no h daughters, sources before/after the pulse, angular collision kernel, interactions, actual core spatial correlation, or derived formation.'}
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
