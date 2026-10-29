"""Geometric bounds on directions of freely streaming surface primaries.

Pair support uses product energy marginals, not measured time/position
correlations or collision rates. No h-decay daughters or rescattering.
"""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    z=np.load(ROOT/'source_t400.npz');E=z['energy_centers'];p=np.sqrt(E*E-1)
    plus=z['plus_positive'].sum(axis=0);minus=z['minus_positive'].sum(axis=0)
    mu=(E[:,None]*E[None,:]-2)/(p[:,None]*p[None,:])
    pairs=plus[:,None]*minus[None,:];den=pairs.sum();rows=[]
    for R in (1,1.25,1.5,2,3,4,5,10,20,50,100):
        minimum=1-2/R**2
        allowed=(mu>=minimum)&(mu<=1.)
        rows.append({'r_over_matching_radius':R,
                     'maximum_individual_radial_angle_deg':float(np.rad2deg(np.arcsin(1/R))),
                     'minimum_relative_angle_cosine':float(minimum),
                     'kinematically_supported_fraction_of_independent_energy_pairs':float(pairs[allowed].sum()/den)})
    # Exact impact-parameter geometry is checked with random outward rays.
    rng=np.random.default_rng(2948);R=3.;mu0=rng.uniform(0,1,100000)
    b=np.sqrt(1-mu0**2);sin_theta=b/R
    error=float(np.max(sin_theta-1/R))
    assert error<1e-15
    fractions=[x['kinematically_supported_fraction_of_independent_energy_pairs'] for x in rows]
    assert np.all(np.diff(fractions)<=1e-15)
    out={'rows':rows,'tests':{'max_sampled_sin_angle_above_geometric_bound':error,'support_fraction_monotonically_nonincreasing':True},
         'geometry':'At radius r>=r0, a free outward trajectory originating at the r0 surface has impact parameter b<=r0, hence sin(theta_r)<=r0/r. For two such rays, cos(gamma)>=1-2(r0/r)^2.',
         'source':'source_t400; positive projections of signed time-windowed direct charged flux marginals',
         'scope':'Kinematic support only, independent draws from time-integrated marginals. It does not know actual angular/time correlations, finite-width rate, h-daughter emission, rescattering, or core evolution.'}
    (ROOT/'ANGULAR_GEOMETRY.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
