"""Same-author controls for the explicit changed-potential candidate."""
import json,math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from match import ROOT,RF,WORK,old


def main():
    read=lambda n:json.loads((ROOT/n).read_text())
    checks=[]
    def group(name,condition,diagnostics):checks.append({'name':name,'pass':bool(condition),'diagnostics':diagnostics})
    point=read('CANDIDATE_M1.json');d=read('DEFORMATION_RESULTS.json');r=read('DEFORMED_MATCHED_RESULTS.json');c=read('CONSEQUENCES.json')
    p=point['heavy'];lam=.03;eps=.002
    u=np.linspace(-20,20,2001)
    potential=lam/4*(u*u-1)**2-eps*(u-1)**3
    factored=(u-1)**2*(lam/4*(u+1)**2-eps*(u-1))
    disc=(lam-3*eps)**2-12*lam*eps
    ly_eff=1-27*p['kappa']**4/(16*p['lambda_X']**3)
    # Lower bound at |Sigma|>=M²/(2mu), after minimizing over X and Y.
    m=30.;mu=old.MU;h=40.;threshold=m*m/(2*mu);ratio2=h*h/(2*threshold**2)
    positivity_margin=(1-ratio2)**2*threshold**2-mu*mu/ly_eff
    pe=float(np.max(abs(potential-factored))/max(1,np.max(potential)))
    group('one_vacuum_and_hidden_tree_bound',disc<0 and 0<eps<2*lam and ly_eff>0 and positivity_margin>0 and pe<1e-14,
          {'stationary_quadratic_discriminant':disc,'factorization_relative_error':pe,
           'remaining_Y_quartic_lower_bound':ly_eff,'large_Sigma_positivity_margin':positivity_margin})
    primary=d['primary'];coarse=next(x for x in d['rows'] if x['tag']=='grid_0.002')
    mass_change=abs(primary['mass_GeV']/coarse['mass_GeV']-1)
    group('new_profile_binding_first_law',primary['classical_absolute_binding'] and primary['Omega']<1
          and abs(primary['charge_relative_residual'])<1e-10 and abs(primary['virial_relative_residual'])<1e-9
          and mass_change<1e-8 and abs(d['first_law']['relative_difference'])<1e-6 and d['first_law']['neighbor_dOmega_dQ']<0,
          {'binding_fraction':1-primary['E_over_Qm'],'domain_tolerance_mass_change':mass_change,
           'neighbor_first_law_relative_error':d['first_law']['relative_difference']})
    stab=read('STABILITY_RESULTS.json')['rows']
    # Known exact symmetries plus expected second-order spatial errors.
    t0=stab[0]['ell1_lowest_translation_expected_zero'][0]
    t1=stab[1]['ell1_lowest_translation_expected_zero'][0]
    z0=stab[0]['phase_lowest_U1_expected_zero'][0];z1=stab[1]['phase_lowest_U1_expected_zero'][0]
    extrap_t=(4*t1-t0)/3;extrap_z=(4*z1-z0)/3
    group('finite_grid_energy_and_symmetry_modes',all(x['radial_fixed_charge_lowest'][0]>0 for x in stab)
          and all(x['ell1_lowest_translation_expected_zero'][1]>0 for x in stab)
          and 3.8<t0/t1<4.2 and 3.8<z0/z1<4.2 and abs(extrap_t)<1e-6 and abs(extrap_z)<1e-6,
          {'minimum_fixed_charge_radial_eigenvalue':min(x['radial_fixed_charge_lowest'][0] for x in stab),
           'translation_zero_second_order_extrapolation':extrap_t,'phase_zero_second_order_extrapolation':extrap_z,
           'certified_full_continuum_or_gauge_spectrum':False})
    cm=read('COLD_MODE_RESULTS.json')['rows'];err=0.
    for x in cm:
        F,S,om,D,k=x['supplied_F'],x['self_consistent_S'],x['Omega'],x['neutral_curvature_D'],x['peak_k_over_mX']
        stiff=np.array([[k*k+.2*F*F,0.,2*F*S],[0.,k*k,0.],[2*F*S,0.,k*k+D]])
        coriolis=np.array([[0.,2*om,0.],[-2*om,0.,0.],[0.,0.,0.]])
        gen=np.block([[np.zeros((3,3)),np.eye(3)],[-stiff,coriolis]])
        growth=float(np.linalg.eigvals(gen).real.max())
        err=max(err,abs(growth-x['peak_growth_over_mX']))
    group('cold_coupled_instability',err<1e-7 and next(x for x in cm if x['supplied_F']==1.)['peak_growth_over_mX']>.1,
          {'polynomial_vs_generator_growth_error':err,'onset_population_derived':False})
    central=next(x for x in r['rows'] if x['tag']=='M1_b3_g0.5_matched')
    compact=np.load(ROOT/'BACKGROUND_M1.npz')['b3_g0.5'];target_times=np.array([0.,.1,.3,.5,1.,2.])
    init=np.array(central['initial_state']);eta2=(central['physical_H_GeV']/old.MRED)**2
    independent=solve_ivp(lambda t,y:old.rhs(t,y,3.,.5,old.XI,eta2),(0,2),init,method='Radau',rtol=2e-10,atol=1e-12,t_eval=target_times)
    assert independent.success
    idx=[int(np.argmin(abs(compact[:,0]-t))) for t in target_times]
    expected=compact[idx,29:];actual=independent.y.T
    ode_error=float(np.max(abs(actual-expected)/(1+abs(expected))))
    max_gauss=max(x['checks']['max_Gauss'] for x in r['rows']);max_f=max(x['checks']['max_Friedmann'] for x in r['rows'])
    group('physical_reheating_and_independent_method',ode_error<1e-8 and max_gauss<1e-10 and max_f<1e-7
          and all(abs(x['target_mass_ratio']-1)<1e-8 for x in r['rows'] if 'target_mass_ratio' in x),
          {'short_Radau_vs_C_scaled_error':ode_error,'max_Gauss':max_gauss,'max_Friedmann':max_f,
           'mass_match_is_calibration':True})
    reference=read('MATCHED_RESULTS.json');base=next(x for x in reference['rows'] if x['tag']=='b3_g0.5_baseline')
    scale_error=max(abs(x['absolute_yield_magnitude']/(base['absolute_yield_magnitude']*x['uniform_scale']**1.5)-1)
                    for x in reference['rows'] if x['tag'].startswith('scale_'))
    tight=next(x for x in r['rows'] if x['tag']=='M1_tighter')
    tol_error=abs(tight['absolute_yield_magnitude']/central['absolute_yield_magnitude']-1)
    phase=r['phase_response'];slope_error=abs(phase['half_step_dlnY_dphase']/phase['dlnY_dphase']-1)
    group('scale_calibration_tolerance_phase',scale_error<1e-7 and tol_error<1e-7 and slope_error<1e-5,
          {'H_power_scaling_relative_error':scale_error,'tight_yield_relative_change':tol_error,
           'phase_derivative_half_step_relative_change':slope_error})
    loss=read('COLD_LOSS_RESULTS.json');coh=loss['coherent'];born=loss['vacuum_Born_pairs']
    coherent_change=abs(coh[-1]['initial_fractional_core_charge_loss_rate_per_second_for_c1']/coh[0]['initial_fractional_core_charge_loss_rate_per_second_for_c1']-1)
    pair_change=abs(born[-1]['pair_fractional_core_rate_per_second']/born[0]['pair_fractional_core_rate_per_second']-1)
    group('changed_profile_included_charge_loss',loss['included_initial_rate_times_13_8Gyr']<.01
          and coherent_change<1e-4 and pair_change<1e-4 and max(x['source_work_relative_residual'] for x in coh)<1e-6,
          {'included_rate_times_age':loss['included_initial_rate_times_13_8Gyr'],
           'coherent_radius_tolerance_change':coherent_change,'Born_quadrature_change':pair_change,
           'complete_quantum_lifetime':False})
    # Explicit forced-neutral control, normalized to factor out the tiny new force.
    early=compact[compact[:,0]<=20.000001];t=early[:,0];Hf=CubicSpline(t,early[:,2]);xr=CubicSpline(t,early[:,31]);xi=CubicSpline(t,early[:,32])
    forced=solve_ivp(lambda tau,y:[y[1],3-3*Hf(tau)*y[1]-.02*(xr(tau)**2+xi(tau)**2)*y[0]],
                     (0,20),[0.,0.],method='DOP853',rtol=2e-10,atol=1e-12,max_step=.02,dense_output=True)
    force_size=phase['early_background_force_epsilon_v3_over_H3']
    forced_phi_over_H=float(np.max(abs(forced.sol(np.linspace(0,20,1001))[0]))*force_size)
    group('early_deformation_and_clock_static_controls',forced.success and forced_phi_over_H<1e-25
          and all(abs(x['Gamma_over_H_uniform_scale_relative_change'])<1e-12 for x in c['clocks'])
          and c['unscreened_profile_energy'][-1]['relative_to_binding']<1e-6,
          {'max_linear_forced_phi_over_H_through_tau20':forced_phi_over_H,
           'gaussian_and_nonlinear_deformation_error_not_certified':True,
           'unscreened_energy_over_binding_e1e6':c['unscreened_profile_energy'][-1]['relative_to_binding']})
    result={'same_author':True,'groups':checks,'all_pass':all(x['pass'] for x in checks)}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'groups':len(checks),'all_pass':result['all_pass'],'failed':[x['name'] for x in checks if not x['pass']]}))
    if not result['all_pass']:raise SystemExit(1)


if __name__=='__main__':main()
