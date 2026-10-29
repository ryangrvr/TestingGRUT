"""Dimensionful scale checks and bounded hypotheses for the M1 synthesis."""
from pathlib import Path
import json,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent

def main():
 c=json.loads((RF/'grut_matched_candidate_20261027'/'CANDIDATE_M1.json').read_text());aud=json.loads((RF/'grut_m1_core_matching_20261029'/'SCENARIO_AUDIT.json').read_text())
 m=c['light_renormalized_targets']['mX_GeV'];T=.020;g=10.75;Mpl=2.435e18
 H=math.sqrt(math.pi**2*g/90)*T*T/Mpl;D=20000.;ratio=H/m;u=ratio*D
 n=aud['ncore_over_mX3'];d=n**(-1/3);Rws=(3/(4*math.pi*n))**(1/3);Rm=60.
 # For independent uniformly distributed, freely streaming particles, this
 # is the swept cylinder per target volume, not a coherent-field theorem.
 swept=math.pi*Rm*Rm*n*D;candidate=.7*swept
 overlap=1-math.exp(-n*4*math.pi/3*(D+2*Rm)**3)
 q=np.load(RF/'grut_m1_open_stress_20261029'/'field_diagnostics.npz');s=q['pulse0025_stress'];t=s['tau'];hist=q['pulse0025_history'];E0=float(hist['E'][0]*100)
 pressures={}
 for r in (20,60):
  tr=s[f'trace{r}_center'];pressures[str(r)]={'signed_integral_mX_time':float(np.trapezoid(tr/3,t)),
     'absolute_integral_mX_time_sampled':float(np.trapezoid(abs(tr/3),t)),
     'absolute_expansion_work_mX_sampled':float(ratio*np.trapezoid(abs(tr),t)),
     'max_abs_pressure_mX_sampled':float(np.max(abs(tr/3)))}
 before=aud['energy_decompositions_at_scale_factor_10']['parent'];after=aud['energy_decompositions_at_scale_factor_10']['through_20000']
 cold_diff=(after['cold_excess_percent']-before['cold_excess_percent'])*aud['mass_reference_mX']/100
 released=before['held_remainder_mX']-after['held_remainder_mX'];f=cold_diff/released
 hypothetical_hot=after['cold_excess_percent']+after['gas_kinetic_percent']+f*after['held_remainder_percent']
 hypothetical_cold=after['cold_excess_percent']+f*after['held_remainder_percent']
 shifts={}
 for label,frac in [('through_20000',after['full_percent']/100),('hypothetical_kept_gas_heat',hypothetical_hot/100),('hypothetical_fully_cooled',hypothetical_cold/100)]:
  scale=(1+frac)**(-2/3);shifts[label]={'mass_increase_fraction':frac,'H_ratio_if_H_three_halves_yield_scaling':scale,
    'H_shift_percent':100*(scale-1),'H_new_GeV':c['heavy']['H_GeV']*scale}
 # Vphi=(s-1)^2*(.75*s^2+1.3*s+.95) and its quadratic discriminant is <0.
 out={'inflation_H_GeV':c['heavy']['H_GeV'],'mX_GeV':m,'T_GeV':T,'gstar':g,'H_at_T_GeV':H,'H_over_mX_at_20MeV':ratio,
  'duration_mXt':D,'delta_ln_a_constant_H':u,'canonical_closed_volume_fractional_abs_work_upper':math.expm1(3*u),
  'initial_total_field_energy_mX':E0,'canonical_energy_gronwall_abs_work_upper_mX':E0*math.expm1(3*u),
  'sampled_regional_flat_pressure_diagnostics':pressures,
  'positivity_proof':{'neutral_potential_factor':'(s-1)^2*(0.75*s^2+1.3*s+0.95)','quadratic_discriminant':1.3**2-4*.75*.95,
      'field_potential':'neutral term + 0.5*s^2*(x^2+y^2)+0.025*(x^2+y^2)^2','claim':'In minimally coupled canonical scalar model V>=0, |trace(Tij)/3|<=K+G+V pointwise. For a comoving volume without boundary flux, |dE/dt|<=3HE and Gronwall gives the stated bound. Curvature improvement, boundary transport, interactions and H evolution require their own treatment.'},
  'density':{'ncore_over_mX3':n,'d_mean_separation_mX_inverse':d,'spherical_wigner_seitz_radius_mX_inverse':Rws,
    'volume_per_core_mX_inverse_cubed':1/n,'probe_radius_mX_inverse':Rm,'maximum_swept_uniform_free_streaming_probability':swept,
    'v_0p7_swept_probability':candidate,'poisson_probability_one_center_within_causal_distance':overlap,
    'scope':'Cylinder estimate assumes uniform independent particle bath, ballistic trajectories, geometric first crossing and v<=1; it does not bound coherent focusing, correlations, re-encounters, long-range scattering or the actual physical mean density at this time.'},
  'calibration':{'baseline_H_GeV':c['heavy']['H_GeV'],'yield_proportional_to_H_power':1.5,'scenarios':shifts,
    'scope':'Sensitivity under fixed dimensionless early history, fixed light mass and entropy factors, with final composition supplied; no formation or physical abundance prediction.'}}
 out['review_arithmetic_correction']={'prior_two_trajectory_cold_per_remainder_ratio':f,
   'hypothetical_if_current_gas_heat_is_retained_percent':hypothetical_hot,
   'hypothetical_if_current_gas_heat_also_fully_redshifts_percent':hypothetical_cold,
   'current_gas_kinetic_percent':after['gas_kinetic_percent'],
   'scope':'Both are hypothetical. The earlier synthesis called 0.841% an arithmetic error without considering complete later cooling; that statement was too broad. Neither value is a dynamical prediction.'}
 (ROOT/'scale_checks.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
