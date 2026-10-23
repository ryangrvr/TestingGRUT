"""Best point-WIMP imitation of the reference reconstructed recoil shape.

This signal-only Kullback--Leibler comparison is an optimistic information
diagnostic. It is not a hypothesis test using a detector's backgrounds.
"""
import json
import numpy as np
from scipy.integrate import simpson
from scipy.optimize import minimize_scalar
from model import ROOT,Object,bin_integrals


def run():
 pub=np.load(ROOT/'data/public_inputs.npz')
 det=json.loads((ROOT/'DETECTOR_RESULTS.json').read_text())
 ref=next(x for x in det['cases'] if x['Q']==39944 and x['v_phi_GeV']==1.)
 obj=Object(39944)
 raw=np.load(ROOT/'data/reconstructed_spectra.npz')['Q39944_v1_future_reconstructed']
 p=raw/sum(raw);mask=p>0
 def point(mass):
  fun=lambda e:obj.rate(e,1.,standard_mass=mass,standard_sigma=1e-45)
  return 20*bin_integrals(fun,pub['future_response_xedges'])@pub['future_response_values']
 def divergence(logm):
  r=point(np.exp(logm));q=r/sum(r)
  return float(np.sum(p[mask]*np.log(p[mask]/np.maximum(q[mask],1e-300))))
 grid=np.linspace(np.log(3.),np.log(1e5),300)
 vals=np.array([divergence(x) for x in grid]);j=int(np.argmin(vals))
 fit=minimize_scalar(divergence,bounds=(grid[max(j-1,0)],grid[min(j+1,len(grid)-1)]),method='bounded',options={'xatol':1e-10})
 mass=float(np.exp(fit.x));r=point(mass);q=r/sum(r)
 count=float(sum(raw)*(ref['kappa_anchor']/1e-9)**2)
 cross=1e-45*count/sum(r)
 result={'scope':'Free mass and normalization point-WIMP imitation in the historical XENONnT 20 ty response. Signal-only ideal information; no background likelihood or claimed significance.',
  'reference_Q':39944,'reference_v_phi_GeV':1.,'reference_kappa':ref['kappa_anchor'],
  'reference_object_mass_GeV':ref['mass_GeV'],'reference_expected_signal_events':count,
  'best_point_WIMP_mass_GeV':mass,'best_point_WIMP_sigma_n_cm2':cross,
  'KL_object_to_best_point_per_signal_event':float(fit.fun),
  'twice_expected_ideal_log_likelihood_ratio_at_reference':float(2*count*fit.fun),
  'signal_events_for_ideal_twice_mean_log_ratio_9':float(9/(2*fit.fun)),
  'fit_success':bool(fit.success),'total_variation_distance':float(.5*np.sum(abs(p-q)))}
 # An alternative LZ approximation: imitate the *soft* efficiency-weighted
 # shape with a free point-WIMP mass, instead of rescaling the hard WIMP
 # limit at the object's actual mass. The disagreement measures a limitation
 # of count rescaling; neither approximation is the full LZ likelihood.
 ee=np.linspace(.001,100,4001)
 eff=np.interp(ee,pub['lz_efficiency_x'],pub['lz_efficiency_y0'],left=0,right=0)
 ob=obj.rate(ee,1.,kappa=ref['kappa_anchor'])*eff
 norm=simpson(ob,x=ee);p_lz=ob/norm
 def lz_point(m):return obj.rate(ee,1.,standard_mass=m,standard_sigma=1e-45)*eff
 def lz_div(logm):
  r=lz_point(np.exp(logm));rnorm=simpson(r,x=ee)
  if rnorm<=0:return np.inf
  q=r/rnorm
  # Hellinger distance stays well-defined at the exact point-WIMP recoil
  # endpoint. The true-spectrum KL would diverge for mismatched supports;
  # this approximate screen has no supplied LZ resolution to remove that.
  return float(.5*simpson((np.sqrt(p_lz)-np.sqrt(q))**2,x=ee))
 lv=np.array([lz_div(x) for x in grid]);lj=int(np.argmin(lv))
 lf=minimize_scalar(lz_div,bounds=(grid[max(lj-1,0)],grid[min(lj+1,len(grid)-1)]),method='bounded')
 lm=float(np.exp(lf.x));ls=1e-45*norm/simpson(lz_point(lm),x=ee)
 llim=float(np.exp(np.interp(np.log(lm),np.log(pub['lz_si_x']),np.log(pub['lz_si_y0']))))
 result['LZ_alternative_shape_matched_ESTIMATE']={
  'best_point_mass_GeV':lm,'apparent_point_sigma_at_anchor_cm2':ls,
  'official_LZ_point_limit_at_fitted_mass_cm2':llim,
  'kappa90_ESTIMATE':float(ref['kappa_anchor']*np.sqrt(llim/ls)),
  'efficiency_weighted_squared_Hellinger_distance':float(lf.fun),
  'ratio_to_same_mass_count_reweighting_ESTIMATE':float(ref['kappa_anchor']*np.sqrt(llim/ls)/ref['lz2024_count_reweighted_kappa90_ESTIMATE'])}
 (ROOT/'SHAPE_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 edges=pub['future_response_yedges'];centers=(edges[1:]+edges[:-1])/2
 np.savetxt(ROOT/'data/reference_shape_comparison.csv',np.column_stack([centers,p,q]),delimiter=',',header='reconstructed_energy_keV,object_probability,best_point_WIMP_probability',comments='')
 print(json.dumps(result,indent=2))


if __name__=='__main__':run()
