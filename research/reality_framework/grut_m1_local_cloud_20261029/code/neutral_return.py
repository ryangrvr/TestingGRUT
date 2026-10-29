"""Ballistic vacuum-decay return of the supplied direct neutral radiation.

This computes a well-defined contribution, not the interacting cloud. Neutral
parents are assigned outward straight radial trajectories and decay according
to the vacuum width. The daughter crossing count is not an absorption or
conversion probability and no rescattering is included.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_expansion_clock_20261028'
MH=np.sqrt(6.)
GAMMA=.00018756589919939715
BETA_STAR=np.sqrt(1-4/MH**2)
E_STAR=MH/2
OMEGA=.8817204919524214


def energy_angle(Eh, x):
    """Charged-daughter lab energy and radial momentum; x=cos(theta*)."""
    gamma=Eh/MH
    beta=np.sqrt(max(0.,1-1/gamma**2))
    return gamma*E_STAR*(1+beta*BETA_STAR*x),gamma*E_STAR*(beta+BETA_STAR*x)


def backward_fraction(Eh):
    gamma=Eh/MH; beta=np.sqrt(max(0.,1-1/gamma**2))
    return max(0.,1-beta/BETA_STAR)/2


def circle_fraction(Eh, radius, emit_radius):
    """Fraction of one daughter's isotropic rest angles hitting inner sphere.

    No classical-deflection approximation is used for the intersection. The
    parent must be radially outward; space is flat during the microscopic event.
    """
    if emit_radius < radius:raise ValueError('Emission point inside target sphere')
    gamma=Eh/MH
    beta=np.sqrt(max(0.,1-1/gamma**2));B=beta/BETA_STAR
    if B>=1:return 0.
    if emit_radius==radius:return (1-B)/2
    D=gamma**2*radius**2/(emit_radius**2-radius**2)
    # Algebraically rationalized to avoid cancellation at small solid angle.
    return .5*D*(1-B)**2/(1+D*(1-B)+np.sqrt(1+D*(1-B**2)))


def return_probability(Eh, radius, emit_radius):
    if Eh<MH:raise ValueError('Parent below vacuum mass shell')
    if emit_radius<radius:raise ValueError('Emission point inside target sphere')
    if Eh>=3.:return 0.
    ell=np.sqrt(max(0.,Eh**2-MH**2))/(MH*GAMMA)
    if ell==0:return circle_fraction(Eh,radius,emit_radius)
    # u is an exponential proper path variate; integrating to infinity includes
    # decays after the archived field simulation. Times are not synchronized.
    f=lambda u:np.exp(-u)*circle_fraction(Eh,radius,emit_radius+ell*u)
    return quad(f,0,np.inf,epsabs=2e-12,epsrel=2e-9,limit=200)[0]


def source(tag):
    out=[];pins={};audit={}
    for window in ('0_1000','1000_10000'):
        p=PARENT/'sources'/f'ringdown_{tag}_spectrum_{window}.npz'
        pins[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        d=np.load(p);counts=d['neutral_N_signed'];energies=d['neutral_EN_signed']
        neg=counts<0;positive=counts>0
        E=np.divide(energies,counts,out=np.zeros_like(counts),where=positive)
        valid=positive&(E>=MH)&np.isfinite(E)
        audit[window]={'signed_number':float(counts.sum()),'discarded_negative_signed_number':float(counts[neg].sum()),
                       'discarded_positive_invalid_number':float(counts[positive&~valid].sum()),
                       'internal_negative_weight':float(d['neutral_N_negative_part'].sum()),
                       'signed_energy_mX':float(energies.sum())}
        out.extend(zip(E[valid],counts[valid]))
    return np.array(out),audit,pins


def checks():
    # At rest the angular integral reduces to the exact subtended solid angle.
    errors=[]
    for r in (20.001,25.,40.,60.,100.,1000.):
        exact=(1-np.sqrt(1-(20/r)**2))/2
        errors.append(abs(circle_fraction(MH,20,r)-exact))
    # Compare against a separate root of the boosted-momentum intersection.
    for E in (2.46,2.6,2.8,2.95,2.999):
        for r in (21.,30.,60.,200.):
            gamma=E/MH;beta=np.sqrt(1-1/gamma**2);B=beta/BETA_STAR
            def cos_lab(x):
                radial=gamma*(beta+BETA_STAR*x)
                transverse=BETA_STAR*np.sqrt(max(0.,1-x*x))
                return radial/np.sqrt(radial**2+transverse**2)
            cutoff=-np.sqrt(1-(20/r)**2)
            x=brentq(lambda x:cos_lab(x)-cutoff,-1.,-B,xtol=1e-14)
            errors.append(abs(circle_fraction(E,20,r)-(x+1)/2))
    # Forward/back daughters obey total energy and momentum in every sample.
    conservation=[];bounds=[]
    for E in (MH,2.55,2.8,2.99,3.,4.):
        g=E/MH;b=np.sqrt(max(0.,1-1/g**2))
        for x in np.linspace(-1,1,201):
            a,pa=energy_angle(E,x);c,pc=energy_angle(E,-x)
            conservation.append(max(abs(a+c-E),abs(pa+pc-np.sqrt(max(E*E-6,0)))))
            if pa<0:bounds.append(a)
    max_backward=max(bounds)
    out={'solid_angle_and_independent_root_max_absolute_error':max(errors),
         'two_daughter_energy_and_radial_momentum_max_absolute_error_mX':max(conservation),
         'largest_sampled_inward_daughter_energy_mX':max_backward,
         'analytic_inward_energy_upper_bound_mX':E_STAR,
         'min_positive_reverse_exchange_threshold_mX':1+2*OMEGA,
         'energies_above_three_can_be_backward':False}
    assert max(errors)<2e-12 and max(conservation)<2e-12 and max_backward<=E_STAR+1e-12
    assert E_STAR < 1+2*OMEGA
    return out


def main():
    verified=checks();rows=[]
    for tag in ('h01','h005','h0025'):
        spec,audit,pins=source(tag);E=spec[:,0];n=spec[:,1]
        backward=np.array([backward_fraction(x) for x in E])
        cases=[]
        for emit in (20.,25.,40.,60.):
            ret=np.array([return_probability(e,20.,emit) for e in E])
            # Per-species equality follows charge-conjugate two-body decay.
            cases.append({'emission_radius_mX_inverse':emit,'target_radius_mX_inverse':20.,
                          'mean_crossings_per_species':float(n@ret),
                          'mean_both_species_crossings':float(2*n@ret),
                          'all_incident_antiparticles_convert_max_direct_exchange_charge_loss':float(2*n@ret)})
        rows.append({'source':tag,'neutral_count':float(n.sum()),'neutral_energy_mX':float(E@n),
                     'neutral_number_below_three_mX':float(n[E<3].sum()),
                     'mean_direct_inward_per_species_any_impact':float(n@backward),
                     'cases':cases,'source_audit':audit,'input_sha256':pins})
    result={'rows':rows,'checks':verified,'width_over_mX':GAMMA,'daughter_rest_speed':BETA_STAR,
            'scope':'Direct vacuum decays of the measured finite-window neutral source; radial parent rays, fixed target sphere, no interactions, absorption, local cloud, or late emission after mXt=10000. Emission at r<60 is a source-radius sensitivity assumption, not measured flux at that radius. Crossing is not capture.'}
    (ROOT/'NEUTRAL_RETURN.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
