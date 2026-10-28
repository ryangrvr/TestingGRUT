"""Carry the same calibrated scale into matching, clocks and future failure tests."""
from pathlib import Path
import json,math
import numpy as np
from scipy.integrate import simpson,cumulative_trapezoid

ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent


def save(name,obj):
    (ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def main():
    old=json.loads((RF/'grut_heavy_vector_budget_20261025/INPUTS.json').read_text())
    data=json.loads((ROOT/'DEFORMED_MATCHED_RESULTS.json').read_text())
    reference=json.loads((ROOT/'MATCHED_RESULTS.json').read_text())
    deformation=json.loads((ROOT/'DEFORMATION_RESULTS.json').read_text())
    row=next(x for x in data['rows'] if x['tag']=='M1_b3_g0.5_matched')
    base=next(x for x in reference['rows'] if x['tag']=='b3_g0.5_baseline')
    s=row['uniform_scale'];p=old['normalized_parameters'];l=old['light_sector'];b=old['benchmark']
    heavy={key:float(p[key]*s) for key in ('H_GeV','M_Y_GeV','h_GeV','mu_GeV')}
    H,M,h,mu=(heavy[key] for key in ('H_GeV','M_Y_GeV','h_GeV','mu_GeV'))
    r=np.sqrt(2)*mu*h/M**2;B=p['kappa']**2/(M*M*(1-r*r));A=r*B/2
    clocks=[]
    for e in (1e-7,1e-6):
        ma=6*e*h;gam=e*e*ma/(48*np.pi)*(1-4*l['mX_GeV']**2/ma**2)**1.5
        old_ma=6*e*p['h_GeV'];old_gam=e*e*old_ma/(48*np.pi)*(1-4*l['mX_GeV']**2/old_ma**2)**1.5
        C=row['radiation_C'];a_mean=np.sqrt(np.pi*C*H/(2*gam))
        Tmean=(90/(np.pi**2*106.75))**.25*np.sqrt(C*2.435e18*H)/a_mean
        clocks.append({'e':e,'mA_GeV':float(ma),'GammaA_GeV':float(gam),
                       'median_tau':float(np.log(2)*H/gam),'mean_decay_scale_factor':float(a_mean),
                       'Tvis_at_mean_decay_GeV_gstar106_75':float(Tmean),
                       'Gamma_over_H_uniform_scale_relative_change':float((gam/H)/(old_gam/p['H_GeV'])-1),
                       'actual_primordial_vector_population_recomputed':False})
    # Conservative unscreened energy for the existing spherical charge profile.
    prof=np.loadtxt(ROOT/'PROFILE_M1.csv',delimiter=',',skiprows=1)
    rho,F,Fp,S,Sp=prof.T;radial=cumulative_trapezoid(rho*rho*F*F,rho,initial=0)
    enclosed=l['Q_reference']*radial/radial[-1]
    integral=simpson(enclosed[1:]**2/rho[1:]**2,x=rho[1:])+l['Q_reference']**2/rho[-1]
    energies=[]
    epq=data['object_energy_per_charge_GeV'];binding=l['Q_reference']*(l['mX_GeV']-epq)
    for e in (1e-7,1e-6):
        energy=e*e*l['gX']*l['vPhi_GeV']*integral/(8*np.pi)
        energies.append({'e':e,'unscreened_E_GeV':float(energy),'relative_to_binding':float(energy/binding)})
    # Leading isolated complex-scalar quartic running: an EFT-range diagnostic only.
    cutoff=1e15;coef=20/(16*np.pi**2);logrange=math.log(cutoff/M)
    quartic=1/(1-coef*logrange)
    reaches4pi=M*math.exp((1-1/(4*np.pi))/coef)
    family=[]
    for eta in (.2,.35,.5,.8,1.):
        req=b['DM_per_visible_entropy_GeV']/(eta*epq+(1-eta)*l['mX_GeV'])
        scale=(req/base['absolute_yield_magnitude'])**(2/3)
        family.append({'supplied_eta':eta,'surviving_pairs_assumed':0.,'required_YQ':req,
                       'scale_prediction_before_new_ODE_run':scale,'H_prediction_GeV':p['H_GeV']*scale,
                       'dynamically_rerun':eta==1.})
    # Fixed future benchmark, not an inferred error bar on a future formation run.
    ratio_m=l['mX_GeV']/epq
    allowance=.01
    test={'frozen_relative_mass_tolerance':allowance,
          'definition':'1+(1-eta)*(mX/epsilonQ-1)+2*(mX/epsilonQ)*p for fixed M1 yield',
          'coefficient_unbound_net_charge':ratio_m-1,'coefficient_pair_ratio':2*ratio_m,
          'minimum_eta_at_p_zero':1-allowance/(ratio_m-1),
          'maximum_p_at_eta_one':allowance/(2*ratio_m),
          'is_observational_confidence_interval':False,'actual_eta':None,'actual_p':None}
    point={'name':'GRUT-M1','status':'specified conditional EFT candidate; single-vacuum binding and repaired mass normalization',
           'parent':'1fcc6b7d37a2c18c8756af5efaed7e30cee461d0',
           'heavy':heavy|{'kappa':p['kappa'],'lambda_X':p['lambda_X'],'lambda_Y':1.,'lambda_Sigma':1.,'e':1e-6},
           'light_renormalized_targets':l,
           'binding_potential_deformation':{'epsilon':.002,'Vphi':'lambda_phi/4*(phi^2-v^2)^2-epsilon*v*(phi-v)^3',
                                           'vacuum_and_curvature_at_plus_v_unchanged':True,
                                           'one_stationary_vacuum':True,'new_supplied_coupling':True},
           'reheating':{'mI_over_H':3.,'GammaI_over_H':.5,'mI_GeV':3*H,'GammaI_GeV':.5*H,'xiX':-5/6},
           'SM_Higgs_curvature':{'xiH':.1,'status':'supplied stabilizing choice motivated by the prior review; no inflation/reheating survival probability computed here'},
           'preparation':{'initial_fields_H_units':row['initial_state'][2:8],
                          'physical_field_velocities':'zero at the selected initial slice',
                          'Sigma_phase_convention':'-pi/2 relative to real X',
                          'state_is_supplied':True,'physical_phase_holder_derived':False,
                          'phase_power_requirement':data['phase_response']},
           'SM_portals_at_matching':{'epsilon_BX':0.,'b_Hphi':0.,'lambda_HX':0.,'status':'no non-gravitational portal in conservative candidate'},
           'vacuum_density_GeV4_input':2.52e-47,
           'EFT_cutoff_GeV_declared':cutoff,
           'higher_operators':'matching inputs assumed negligible for this conditional calculation; no UV completion derived',
           'other_phi_odd_scalar_operators':'renormalized matching choices at the reference scale, not all-order zeros after dropping phi parity',
           'final_population_matching_hypothesis':{'Q':200000,'eta':1.,'p':0.,'YQ_target':data['target_YQ_for_all_reference_objects_no_pairs'],
                                                  'actual_formation_computed':False},
           'calibrated_quantity':'early physical scale, using the target mass budget',
           'physical_final_abundance':None}
    save('CANDIDATE_M1.json',point)
    save('CONSEQUENCES.json',{'uniform_scale':s,'r':float(r),'A_GeV_minus2':float(A),'B_GeV_minus2':float(B),
         'relative_A_change':s**-2,'relative_included_cold_leakage_rate_change_at_fixed_light_profile':s**-4,
         'clocks':clocks,'unscreened_profile_energy':energies,
         'finite_cutoff_diagnostic':{'cutoff_GeV':cutoff,'isolated_one_loop_complex_quartic_at_cutoff':quartic,
                                    'scale_at_quartic_4pi_GeV':reaches4pi,'full_coupled_RG_completed':False},
         'scalar_source_invariance':'dimensionless regulated scalar approximation only; no full quantum/gauge covariance transfer certified',
         'heat_source_invariance':'dimensionless decay width/H clocks checked; populations and physical freeze-out not recomputed',
         'zero_portal_direct_detection_rate':0.,
         'portal_extension_conditional_node_keV':deformation['primary']['ideal_node_keV'],
         'portal_extension_is_cosmologically_certified':False,'population_family':family,'formation_failure_test':test,
         'scalar_cubic_relative_to_original':(.03-.002)/.03,
         'collision_kernels_recomputed_for_deformation':False})


if __name__=='__main__':main()
