"""Scoped formation supply, two-body seed, condensate onset and wall calculations.

No nonlinear fragmentation, radiation spectrum or object mass function is
inferred from linear growth. All charge and condensate preparation is input.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad,simpson
from scipy.optimize import brentq,minimize_scalar
from scipy.special import zeta

ROOT=Path(__file__).resolve().parents[1]
EMP=ROOT.parent/'grut_portal_empirical_test_20261022'
G,LAM,LX,HD=.1,.01,.001,.01
MPL,HBAR,HBARC=2.435e18,6.582119569e-25,.1973269804
DM_ENTROPY,TI,GSTAR=4.4e-10,.020,10.75
ATTR=G**4/LAM-LX


def condensate(number,v):
 mass=G*v
 f2=number/mass
 for _ in range(12):f2=number/np.sqrt(mass*mass-ATTR*f2)
 omega=np.sqrt(mass*mass-ATTR*f2)
 phi2=v*v-G*G*f2/LAM
 assert phi2>0
 return f2,omega,phi2


def growth(number,v,y):
 """Positive real growth root of the full three-field dispersion.

 y = k^2/(ATTR*f^2). Solve for Gamma^2 normalized by the NR peak,
 avoiding subtraction of nearby O(m_X^2) frequency eigenvalues.
 """
 f2,w,phi2=condensate(number,v)
 a=2*LX*f2;b=2*LAM*phi2;c2=4*G**4*phi2*f2
 k2=ATTR*f2*y
 scale=(ATTR*f2/(2*w))**2
 def det(u):
  zz=scale*u
  return (((k2+a+zz)*(k2+zz)+4*w*w*zz)*(k2+b+zz)-c2*(k2+zz))/(ATTR*ATTR*f2*f2*b)
 if y<=0 or det(0)>=0:return 0.
 hi=2.
 while det(hi)<0:hi*=2
 return float(np.sqrt(scale*brentq(det,0,hi,xtol=1e-14,rtol=1e-14)))


def linear_case(number,v):
 f2,w,phi2=condensate(number,v)
 fit=minimize_scalar(lambda y:-growth(number,v,y),bounds=(.001,1.999),method='bounded',options={'xatol':1e-11})
 peak=growth(number,v,fit.x);k=np.sqrt(ATTR*f2*fit.x)
 mass=G*v;ms=np.sqrt(2*LAM)*v
 # Exact stable formula for the positive k^2 band edge.
 a=2*LX*f2;b=2*LAM*phi2;c2=4*G**4*phi2*f2
 kcrit2=2*(c2-a*b)/(np.sqrt((a-b)**2+4*c2)+a+b)
 wavelength=2*np.pi/k
 hubble=np.sqrt(np.pi**2*GSTAR/90)*TI*TI/MPL
 cell_charge=number*wavelength**3
 cs2=-ATTR*f2/(2*mass*mass-3*ATTR*f2)
 return {'supplied_charge_density_GeV3':number,'f_squared_GeV2':f2,'omega_GeV':w,'phi_squared_GeV2':phi2,
  'cs_squared_from_charge_compressibility':cs2,'maximum_growth_GeV':peak,'maximum_growth_over_H':peak/hubble,
  'growth_time_seconds':HBAR/peak,'fastest_y':float(fit.x),'fastest_k_GeV':k,
  'band_edge_k_GeV':float(np.sqrt(kcrit2)),'wavelength_GeV_inverse':wavelength,
  'wavelength_cm':wavelength*HBARC*1e-13,'charge_in_wavelength_cube_DIAGNOSTIC':cell_charge,
  'nonrelativistic_growth_relative_difference':float(peak/(ATTR*f2/(2*w))-1),
  'ideal_gas_BEC_temperature_GeV_DIAGNOSTIC':float(2*np.pi/mass*(number/zeta(1.5))**(2/3)),
  'coherent_state_shot_noise_to_order_one_time_over_Hubble_DIAGNOSTIC':float(.5*np.log(cell_charge)*hubble/peak),
  'number_conservation_relative_residual':float(w*f2/number-1),'mediator_to_growth_scale_ratio':float(ms/peak)}


def run():
 det=json.loads((EMP/'DETECTOR_RESULTS.json').read_text())
 hist=json.loads((EMP/'HISTORY_RESULTS.json').read_text())
 entropy=2*np.pi**2*GSTAR*TI**3/45
 rows=[]
 for d in det['cases']:
  q,v=d['Q'],d['v_phi_GeV'];mass=d['mass_GeV']
  h=next(r for r in hist['cases'] if r['Q']==q and r['v_phi_GeV']==v)
  yq=DM_ENTROPY*q/mass;number=yq*entropy
  release=G*v*q-mass;budget=release*DM_ENTROPY/mass
  sigma=2*np.sqrt(2*LAM)*v**3/3
  alpha=G**4*v*v/(4*np.pi*(G*v)**2);reduced=G*v/2;ms=np.sqrt(2*LAM)*v
  trace=2*reduced*alpha/ms
  lc=linear_case(number,v)
  row={'Q':q,'v_phi_GeV':v,'supplied_portal_kappa_from_parent':d['kappa_anchor'],
   'reference_object_mass_GeV':mass,'supplied_charge_yield':yq,
   'supplied_object_yield':DM_ENTROPY/mass,
   'empty_cold_portal_max_DM_fraction_nominal':h['cold']['injected_energy_over_reference_DM_energy_BOUND'],
   'empty_cold_portal_max_DM_fraction_enlarged_response':h['cold_production_stress_envelope']['injected_energy_over_reference_DM_energy_BOUND'],
   'free_carrier_mass_GeV':G*v,'mediator_mass_GeV':ms,'NR_Yukawa_alpha':alpha,
   'NR_Yukawa_zero_energy_Birman_Schwinger_trace':trace,
   'two_body_NR_Yukawa_bound_seed_absent':bool(trace<1),
   'free_to_object_release_GeV_per_object':release,
   'binding_release_GeV_per_visible_entropy':budget,'binding_release_over_object_DM_energy':budget/DM_ENTROPY,
   'binding_release_over_prior_portal_budget_nominal':budget/h['cold']['all_energy_survives_matter_GeV_per_entropy_BOUND'],
   'binding_release_over_prior_portal_budget_enlarged':budget/h['cold_production_stress_envelope']['all_energy_survives_matter_GeV_per_entropy_BOUND'],
   'maximum_yield_if_binding_release_as_nonrelativistic_s':budget/ms,
   'radiation_ratio_at_nominal_s_decay_if_all_release_as_s_rest_mass_ESTIMATE':4*budget/(3*h['cold']['decay_temperature_radiation_estimate_GeV']),
   'neutral_wall_tension_GeV3':sigma,'neutral_wall_tension_cuberoot_MeV':1000*sigma**(1/3),
   'wall_tension_over_1MeV_cubed':sigma/1e-9,
   'wall_tension_over_2026_081MeV_cubed':sigma/(.00081**3),'linear_condensate':lc}
  rows.append(row)
 ref=next(r for r in rows if r['Q']==39944 and r['v_phi_GeV']==1.)
 n=ref['linear_condensate']['supplied_charge_density_GeV3']
 ys=np.linspace(0,2.2,441);exact=np.array([growth(n,1.,y) for y in ys])
 f2,w,_=condensate(n,1.);nr=ATTR*f2/(2*w)*np.sqrt(np.maximum(ys*(2-ys),0))
 (ROOT/'data').mkdir(exist_ok=True)
 np.savetxt(ROOT/'data/reference_dispersion.csv',np.column_stack([ys,exact,nr]),delimiter=',',header='k_squared_over_attr_f_squared,full_growth_GeV,NR_growth_GeV',comments='')
 # Analytic kink in u=z/delta: retain dimensionless coordinates for reuse.
 u=np.linspace(-12,12,1201);phi=np.tanh(u);deriv=1/np.cosh(u)**2
 np.savetxt(ROOT/'data/neutral_wall_profile.csv',np.column_stack([u,phi,deriv,deriv**2]),delimiter=',',header='z_over_delta,phi_over_v,dphi_over_v_du,sech4',comments='')
 # Tiny density controls for the full determinant versus compressibility.
 ksmall=1e-5
 full_cs=-(growth(n,1.,ksmall)**2)/(ATTR*f2*ksmall)
 # Kink tension by a distinct direct energy quadrature.
 wall_integral=quad(lambda u:1/np.cosh(u)**4,-20,20,epsabs=1e-12)[0]*LAM*np.sqrt(2/LAM)/2
 # The Yukawa trace integral is separately integrated in x=m_s*r.
 trace_quad=2*(G/2)*(G*G/(4*np.pi))/np.sqrt(2*LAM)*quad(lambda x:np.exp(-x),0,np.inf)[0]
 result={'scope':'Conditional formation prerequisites and linear onset; no final mass function, formation-radiation spectrum or multi-target forecast.',
  'reference_parameters':{'g':G,'lambda_phi':LAM,'lambda_X':LX,'h_D':HD,'T_SM_initial_GeV':TI,'g_star':GSTAR,
    'DM_energy_per_visible_entropy_GeV':DM_ENTROPY,'effective_attractive_quartic':ATTR},
  'cases':rows,'reference_Q':39944,'reference_v_phi_GeV':1.,
  'bulk_Q_matter':{'omega_floor_GeV_at_v1':float((LAM*LX)**.25),'f_squared_over_v_squared':float(np.sqrt(LAM/LX)),
    'phi_zero_branch_consistent':bool(np.sqrt(LAM/LX)>LAM/G**2)},
  'wall_input':{'network_bound_cuberoot_MeV':.81,'source':'Caloni et al., arXiv:2602.20050v1, section 4.2; 95% credible, specified stable-network model',
    'v_max_GeV_for_same_lambda_and_standard_stable_network':float(.00081/(2*np.sqrt(2*LAM)/3)**(1/3)),
    'portal_generates_Z2_bias':False},
  'controls':{'full_small_k_cs_over_thermodynamic_cs_minus_one':float(full_cs/ref['linear_condensate']['cs_squared_from_charge_compressibility']-1),
    'small_k_y_for_control':ksmall,'kink_energy_quadrature_relative_residual':float(wall_integral/ref['neutral_wall_tension_GeV3']-1),
    'Yukawa_trace_quadrature_relative_residual':float(trace_quad/ref['NR_Yukawa_zero_energy_Birman_Schwinger_trace']-1)},
  'literature_and_prior_work':'Methods and source credits in SOURCES.md; older charged-phase/wall work is not rerun. No new primordial charge source is included.'}
 (ROOT/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'reference':ref,'controls':result['controls'],'bulk':result['bulk_Q_matter']},indent=2))


if __name__=='__main__':run()
