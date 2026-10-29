"""Kinematic support for resonant pairs in hypothetical outward cones.

This is not a rate: it uses time-integrated charged flux energy marginals,
without their joint time/position/angular distributions or h-decay daughters.
"""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
WIDTH=.00018756589919939715
MH=np.sqrt(6.)

def source(species,windows):
    n=en=None
    for w in windows:
        d=np.load(ROOT/'sources'/f'ringdown_h0025_spectrum_{w}.npz')
        if n is None:n=d[species+'_N_signed'].copy();en=d[species+'_EN_signed'].copy()
        else:n+=d[species+'_N_signed'];en+=d[species+'_EN_signed']
    negative=float(np.maximum(-n,0).sum());ok=n>0
    energy=np.divide(en,n,out=np.zeros_like(n),where=ok)
    ok &= energy>1
    return energy[ok],n[ok],negative

def run(windows):
    e1,w1,ng1=source('plus',windows);e2,w2,ng2=source('minus',windows)
    p1=np.sqrt(e1*e1-1)[:,None];p2=np.sqrt(e2*e2-1)[None,:]
    e1=e1[:,None];e2=e2[None,:];weight=w1[:,None]*w2[None,:]
    mu=(e1*e2-2)/(p1*p2);delta=MH*WIDTH/(2*p1*p2)
    bands=[]
    for deg in (0,1,5,10,20,30,45,90):
        lower=np.cos(2*np.deg2rad(deg));inside=(mu>=lower)&(mu<=1.)
        within_width=(mu+delta>=lower)&(mu-delta<=1.)
        bands.append({'cone_half_angle_deg':deg,
                      'fraction_of_pairs_for_which_exact_pole_is_geometrically_accessible':float(weight[inside].sum()/weight.sum()),
                      'fraction_of_pairs_that_can_reach_one_width_band':float(weight[within_width].sum()/weight.sum())})
    return {'windows':windows,'positive_number':float(w1.sum()),'negative_number':float(w2.sum()),
            'discarded_negative_binned_weight':{'plus':ng1,'minus':ng2},'bands':bands}

def main():
    # Independent monochromatic geometry: p1=p2=m implies s=6 at 90 degrees;
    # a 20 degree half-cone only permits relative angles up to 40 degrees.
    p=1.;E=np.sqrt(1+p*p);mu=(E*E-2)/(p*p)
    assert abs(mu)<1e-14 and mu<np.cos(np.deg2rad(40))
    results=[run(['0_1000']),run(['0_1000','1000_10000'])]
    for row in results:
        fractions=[x['fraction_of_pairs_for_which_exact_pole_is_geometrically_accessible'] for x in row['bands']]
        assert all(0<=x<=1 for x in fractions) and np.all(np.diff(fractions)>=-1e-15)
    out={'model':'hypothetical common outward cone; both charged momenta lie within half-angle alpha of the same direction',
         'resonant_condition':'mu_res=(Eplus Eminus-2)/(pplus pminus), s=6',
         'collinear_formula':'s=2+2 cosh(yplus-yminus); s=6 requires |delta y|=arcosh(2)',
         'checks':{'monochromatic_mu_res':float(mu),'cone_support_fraction_monotonic':True},
         'results':results,
         'scope':'These are support fractions of time-integrated direct charged flux energy marginals. They are NOT collision rates or source angle measurements. No spatial overlap, scattering angular weight, h-decay daughters, core returns, coherent amplitudes or finite-time spectral uncertainty is modeled.'}
    (ROOT/'BEAM_KINEMATICS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
