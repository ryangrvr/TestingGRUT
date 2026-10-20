"""One fixed cosmological preparation; a falsifiable abundance FLOOR.

python empirical.py predict     # does not open heldout_data.json
python empirical.py compare     # creates prediction, then compares held-out target

The exact net charge, not an invented soliton formation efficiency or lifetime,
closes the conservative exclusion. A predicted soliton mass function remains open.
"""
from pathlib import Path
import hashlib, json, math, sys
from scipy.special import zeta

P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def predict():
    d=json.loads((P/'calibration_data.json').read_text())
    model=json.loads((P/'protocol.json').read_text())['empirical_point']
    g=model['gX'];v=model['v_GeV'];m=model['mX_GeV'];mu=(model['a']*model['b'])**.25
    hbar=d['h_SI']/(2*math.pi);GeV_J=d['eV_J']*1e9;GeV_kg=GeV_J/d['c_SI']**2
    hbar_GeVs=hbar/GeV_J;hbarc_GeVm=hbar_GeVs*d['c_SI']
    MP_kg=math.sqrt(hbar*d['c_SI']/(8*math.pi*d['G_SI']));MP=MP_kg/GeV_kg
    H100=100e3/(1e6*d['parsec_m'])
    rho_c=3*H100**2/(8*math.pi*d['G_SI'])/GeV_kg/1e6 # GeV/cm^3
    n_gamma=2*float(zeta(3))/math.pi**2*(d['kB_SI']*d['T_gamma_K']/(hbar*d['c_SI']))**3/1e6
    # Ideal standard photon+decoupled neutrino entropy, useful only for Ti.
    s0=2*math.pi**2/45*(43/11)*(d['kB_SI']*d['T_gamma_K']/(hbar*d['c_SI']))**3/1e6
    mb=d['mean_mass_per_baryon_GeV_convention']
    nB=d['omega_b']*rho_c/mb;Y=nB/s0
    omega_floor=mu*m*Y*s0/rho_c
    # The same result by a route independent of T_gamma and entropy conversion.
    direct=mu*m/mb*d['omega_b']
    F=.5;S=math.sqrt(.75);Om=math.sqrt(.775);nbar=Om*F*F
    ebar=.5*Om*Om*F*F+.25*(S*S-1)**2+.5*S*S*F*F+.025*F**4
    nQ=g*v**3*nbar;Ti=(nQ/(Y*(2*math.pi**2/45)*d['visible_gs_initial']))**(1/3)
    H=math.sqrt(math.pi**2*d['visible_gstar_initial']/90)*Ti**2/MP
    rhoh=g*g*v**4*ebar;rhor=math.pi**2/30*d['visible_gstar_initial']*Ti**4
    thermal_x_at_T=math.exp(-m/Ti) # displayed to show NOT thermal Boltzmann production
    out={'point':model,'calibration_sha256':sha(P/'calibration_data.json'),'protocol_sha256':sha(P/'protocol.json'),'source_sha256':sha(Path(__file__)),'mass_floor_per_charge_GeV':mu*m,'mu_star':mu,'omega_hidden_floor':omega_floor,'omega_floor_BBN_sigma_only':omega_floor*d['omega_b_sigma']/d['omega_b'],'omega_floor_1percent_mass_systematic':omega_floor*d['mean_mass_systematic_fraction'],'independent_ratio_route':direct,'charge_yield_YQ_equals_late_YB':Y,'present_charge_density_cm3':Y*s0,'rho_critical_over_h2_GeV_cm3':rho_c,'entropy_density_today_cm3':s0,'photon_density_today_cm3':n_gamma,'eta_B':nB/n_gamma,'initial':{'F':F,'S':S,'Omega_bar':Om,'nbar':nbar,'ebar':ebar,'ebar_over_nbar':ebar/nbar,'nQ_GeV3':nQ,'visible_temperature_GeV':Ti,'H_GeV':H,'H_over_mX':H/m,'H_delta_t_for_tau240':240*H/m,'rho_hidden_over_radiation':rhoh/rhor,'H_times_simulated_physical_box':H*64/m,'thermal_boltzmann_factor_NOT_used':thermal_x_at_T},'physical_reference_object':{'mass_GeV':v/g*363.9229985516867,'Q':399.4384624223012/g**2,'R90_m':5.103432213468092/m*hbarc_GeVm,'tau160_seconds':160/m*hbar_GeVs,'mass_kg':v/g*363.9229985516867*GeV_kg},'survival_statement':'Q(a) a^3 conserved for exact adopted hidden U(1), no later visible entropy injection. E>=mu_star*mX*abs(Q) at leading classical/semiclassical weak-gravity order. Fission, evaporation into charged X, and unknown formation efficiency cannot erase this net-charge energy floor. This is NOT a computed individual quantum soliton lifetime.','cosmology_assumptions':['Primordial rotating condensate is an explicitly imposed initial state; no generation of net hidden charge from zero is derived.','YQ is fixed equal to the observed late baryon yield; equality is an added hypothesis, not forced by the action. SM baryogenesis and sphaleron history are not solved.','No hidden charge breaking, late dilution, gravitational collapse hiding global charge, or compensating negative-energy sector.','Leading scalar potential and weak gravitational binding; all-orders quantum/UV completion not established.','For comparison with omega_c, surviving charge carriers form a cold, effectively collisionless late matter component. A full Boltzmann/structure evolution and CMB angular-spectrum likelihood have not been supplied.'],'claim_class':'Conditional conservative abundance constraint for one fully specified initial-yield point; not an unconditional prediction of GRUT or an experimentally verified composite population.'}
    out['checks']={'mass_scale':abs(g*v/m-1)<1e-14,'two_abundance_routes':abs(direct/omega_floor-1)<1e-12,'local_expansion_negligible':240*H/m<1e-5,'radiation_dominated_preparation':rhoh/rhor<1e-6}
    return out
def compare(pred):
    target=json.loads((P/'heldout_data.json').read_text());d=json.loads((P/'calibration_data.json').read_text())
    # No Gaussian significance is quoted: a floor, inferred target and theory
    # assumptions are not a complete likelihood. Use conservative 3-sigma edges.
    conservative_floor=pred['omega_hidden_floor']*(1-3*d['omega_b_sigma']/d['omega_b'])/(1+d['mean_mass_systematic_fraction'])
    upper=target['mean']+3*target['sigma']
    return {'target':target,'prediction_sha256':sha(P/'EMPIRICAL_PREDICTION.json'),'nominal_minimum_overproduction_factor':pred['omega_hidden_floor']/target['mean'],'BBN_3sigma_and_mass_allowance_floor':conservative_floor,'target_3sigma_upper':upper,'conservative_overproduction_factor':conservative_floor/upper,'excluded_under_stated_cosmology':conservative_floor>upper,'derived_constraint_mX_times_YQ_over_YB_GeV':target['mean']*d['mean_mass_per_baryon_GeV_convention']/(pred['mu_star']*d['omega_b']),'no_refit':'The 100 GeV / YQ=YB point is retained and rejected; derived constraint is not fed back to select a new mass.','missing_for_positive_empirical_success':['Derived/calibrated hidden charge production or independently constrained initial yield','Late relic distribution, equation of state and survival under cosmological evolution','Direct forward CMB or another observable likelihood for an independently fixed viable point']}
if __name__=='__main__':
    a=predict();(P/'EMPIRICAL_PREDICTION.json').write_text(json.dumps(a,indent=2)+'\n')
    if not all(a['checks'].values()):raise RuntimeError('compatibility checks failed')
    if len(sys.argv)>1 and sys.argv[1]=='compare':
      b=compare(a);(P/'EMPIRICAL_COMPARISON.json').write_text(json.dumps(b,indent=2)+'\n');print(json.dumps(b,indent=2))
    else:print(json.dumps(a,indent=2))
