#!/usr/bin/env python3
"""A screened-gauge binding-energy bound using admissible charge-neutral data."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import simpson

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
sys.path.insert(0,str(ROOT/'code'))
from pair_rate import transform
from scattering import Background,V,MX


def main():
    bg=Background();m=json.loads((BASE/'grut_z6_gauge_bridge_20261025/MATCHING_RESULTS.json').read_text())
    p=m['parameters'];h=p['Higgs_amplitude_over_M']*p['M_GeV'];e=p['gauge_coupling']
    mass=6*e*h;omega=bg.omega*MX
    book=json.loads((BASE/'grut_portal_debris_window_20261023/FIXED_YIELD_RESULTS.json').read_text())
    row=next(a for a in book['rows'] if a['Q']==200000 and a['lambda_phi']==.03)
    control=[]
    for nr in [6001,12001,18001]:
        r=np.linspace(0,70,nr);F,_=bg.fields(r);Fp=bg.f(r,1)
        q=4*np.pi*omega*V**2/MX**3*simpson(r*r*F*F,x=r)
        n2=4*np.pi*omega**2*V**4/MX**3*simpson(r*r*F**4,x=r)
        grad2=4*np.pi*omega**2*V**4/MX*simpson(r*r*(2*F*Fp)**2,x=r)
        bound=n2/(72*h*h)
        control.append(dict(radial_points=nr,charge=float(q),charge_relative_error=float(q/200000-1),
            integral_density_squared_GeV3=float(n2),gradient_to_density_squared_GeV2=float(grad2/n2),
            added_energy_bound_GeV=float(bound),relative_large_mass_error_bound=float(grad2/(mass*mass*n2))))
    r=np.linspace(0,70,12001);F,_=bg.fields(r)
    q,w=np.polynomial.legendre.leggauss(768);q=(q+1)*10;w=w*10
    ft=transform(r,F*F,q)
    real_norm=4*np.pi*simpson(r*r*F**4,x=r)
    spectral=q*q*ft*ft/(2*np.pi**2)
    parsed=np.sum(w*spectral)
    ratios=[]
    for mu in [.1,1.,10.,1000.,mass/MX]:
        relative=np.sum(w*spectral/(1+(q/mu)**2))/real_norm
        ratios.append(dict(vector_mass_over_mX=mu,screened_energy_to_local_bound=float(relative)))
    best=control[-1];E=row['mass_GeV'];margin=200000*MX-E
    variant=json.loads((ROOT/'POINT_RESULTS.json').read_text())['wall_bias_variant']
    _,S=bg.fields(r)
    bias_energy=4*np.pi/MX**3*variant['DeltaV_GeV4']/4*simpson(r*r*(S-1)**2*(S+2),x=r)
    out=dict(scope='Admissible Gauss-neutral initial data and a positive screened-energy variational bound. Not a complete coupled soliton/Floquet spectrum.',
        controls=control,screened_Fourier_controls=ratios,Parseval_relative_error=float(parsed/real_norm-1),
        matched_free_mass_GeV=MX,old_object_energy_GeV=E,free_particle_binding_margin_GeV=margin,
        added_gauge_energy_bound_GeV=best['added_energy_bound_GeV'],
        bound_over_existing_binding_margin=best['added_energy_bound_GeV']/margin,
        trial_energy_below_free_particles=E+best['added_energy_bound_GeV']<200000*MX,
        proposed_bias_variant_added_trial_energy_GeV=float(bias_energy),
        proposed_bias_plus_gauge_energy_over_binding_margin=float((bias_energy+best['added_energy_bound_GeV'])/margin),
        static_compensating_Higgs_charge=-200000/6,full_current_identity='Q_X+3Q_Y+6Q_Sigma=0; the displayed trial has Y=0.',
        residual_discrete_charge_of_one_reference_core=200000%6,
        minimal_number_of_reference_cores_in_a_neutral_group=3,
        neutral_group_X_charge=600000,neutral_group_Sigma_number_charge=-100000,
        alternative_two_antiX_global_compensation_energy_GeV=2*MX,
        quantum_scope='A single Q=200000 core is 2 mod 6. Classical continuous screening is not a full discrete-gauge quantum-state construction. A compact whole must include compensating sectors; three such cores have total charge divisible by six.',
        minimum_heavy_mass_GeV=min(m['matching']['gauge_vector_mass_GeV'],m['matching']['Higgs_radial_mass_GeV'],*np.sqrt(m['matching']['mediator_squared_masses_GeV2'])),
        full_fluctuation_stability_certified=False,
        assumptions='Low-energy renormalized light parameters matched as in 7737487. The Higgs modulus is h; infinitesimal charge-breaking gives metastability rather than exact fixed-Q stationary states.')
    (ROOT/'SCREENING_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('SCREENING_ENERGY_BOUND_GeV',best['added_energy_bound_GeV'],'BINDING_RATIO',out['bound_over_existing_binding_margin'],'PARSEVAL',out['Parseval_relative_error'])


if __name__=='__main__':main()
