"""Friedmann lower bound: no assumption that hidden matter is cold dust.

The scalar energy-charge floor remains valid on an instantaneous flat FRW
slice: gradient energies scale with a^-2 and stay nonnegative. Conserved
comoving charge and visible entropy fix nQ today. Positive other densities
can only increase H0. No late equation of state or soliton fraction is fitted.
"""
from pathlib import Path
import hashlib, json, math
from empirical import predict

P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prediction():
    d=json.loads((P/'calibration_data.json').read_text());base=predict()
    wb=d['omega_b'];wh=base['omega_hidden_floor'];frac=1-3*d['omega_b_sigma']/wb
    whlow=wh*frac/(1+d['mean_mass_systematic_fraction']);wblow=wb*frac
    return {'point':base['point'],'added_point_specification':'Flat FRW, weak gravitational binding; rho_Lambda=0 (a nonnegative cosmological constant only strengthens the lower bound), nonnegative other densities. Same scalar couplings, scale, initial condensate and yield.','calibration_sha256':sha(P/'calibration_data.json'),'amendment_sha256':sha(P/'EMPIRICAL_AMENDMENT.json'),'source_sha256':sha(Path(__file__)),'formula':'H0 >= 100 sqrt(omega_b+omega_hidden_floor) km/s/Mpc','H0_lower_bound_km_s_Mpc':100*math.sqrt(wb+wh),'H0_lower_bound_BBN_3sigma_mass_allowance':100*math.sqrt(wblow+whlow),'omega_hidden_floor':wh,'charge_yield':base['charge_yield_YQ_equals_late_YB'],'chain':['Primordial rotating nonthermal charged condensate, 3D fragmentation studied at common scalar point; initial hidden asymmetry is an explicit postulate.','YQ=YB(BBN)=8.6758e-11 with conserved visible entropy; no fitted formation efficiency.','Exact adopted U(1) preserves net comoving charge even when particular composites evaporate or split.','Leading scalar energy bound rho_hidden>=mu_star*mX*abs(nQ) on an FRW slice.','Flat-FRW Friedmann constraint implies H0 lower bound independent of hidden pressure.'],'scope':'Conditional exclusion of this preparation at leading classical/semiclassical EFT order, not a positive detection, quantum-lifetime calculation, or an unconditional prediction from the action alone.'}
if __name__=='__main__':
    pred=prediction();pp=P/'GEOMETRIC_PREDICTION.json';pp.write_text(json.dumps(pred,indent=2)+'\n')
    # Prediction is fully constructed before opening comparison data.
    data=json.loads((P/'heldout_megamaser.json').read_text());upper=data['mean']+3*data['sigma']
    d=json.loads((P/'calibration_data.json').read_text());mu=pred['point']['a']**.25*pred['point']['b']**.25
    out={'target':data,'prediction_sha256':sha(pp),'H0_predicted_lower_bound':pred['H0_lower_bound_km_s_Mpc'],'H0_conservative_lower_bound':pred['H0_lower_bound_BBN_3sigma_mass_allowance'],'observed_3sigma_upper':upper,'incompatible_at_conservative_edges':pred['H0_lower_bound_BBN_3sigma_mass_allowance']>upper,'derived_nominal_constraint_mX_times_YQ_over_YB_GeV':((data['mean']/100)**2-d['omega_b'])*d['mean_mass_per_baryon_GeV_convention']/(mu*d['omega_b']),'no_refit':'Keep and reject 100 GeV, YQ=YB, flat nonnegative-energy preparation. Positive Lambda increases predicted H0 minimum.','statistical_limit':'Comparison to published inferred H0 with conservative edges. No Gaussian sigma detection/exclusion significance and no new maser-likelihood analysis claimed.'}
    (P/'GEOMETRIC_COMPARISON.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
