"""Use unmodified ACROPOLIS on an overestimate of the produced scalar population."""
import argparse,json,time,hashlib
import numpy as np
from acropolis.models import DecayModel
import acropolis.flags as flags
import acropolis.params as params
import acropolis.cascade as cascade
import acropolis.nucl as nucl
from model import ROOT,VH
from history import ME,HBAR


def run(only=None):
 flags.verbose=False
 hist=json.loads((ROOT/'HISTORY_RESULTS.json').read_text())
 ref=next(x for x in hist['cases'] if x['Q']==39944 and x['v_phi_GeV']==1.)
 nominal=ref['cold'];stress=ref['cold_production_stress_envelope']
 mass=nominal['s_mass_GeV'];mix=nominal['theta']
 gee=mix*mix*ME*ME*mass/(8*np.pi*VH*VH)*(1-4*ME*ME/(mass*mass))**1.5
 max_tau=HBAR/gee
 cases=[]
 for name,c in [('nominal_all_injected_energy_as_scalars',nominal),('stress_all_injected_energy_as_scalars',stress)]:
  bree=gee/c['s_total_width_in_stated_SM_response_GeV']
  y=c['all_energy_survives_matter_GeV_per_entropy_BOUND']/mass
  cases.append((name,c['lifetime_seconds'],y,bree,1-bree,120,20))
 cases.append(('envelope_stress_yield_photons_at_longest_lifetime',max_tau,
  stress['all_energy_survives_matter_GeV_per_entropy_BOUND']/mass,0.,1.,120,20))
 cases.append(('same_envelope_grid_refinement',max_tau,
  stress['all_energy_survives_matter_GeV_per_entropy_BOUND']/mass,0.,1.,150,30))
 results={'ACROPOLIS_version':'1.3.1','scope':'Incremental photodisintegration on the supplied standard BBN/cosmological background; overestimated scalar populations. Does not refit standard BBN or solve formation.',
          'reference_Q':39944,'reference_v_phi_GeV':1.,'reference_kappa':ref['kappa_anchor'],'cases':[]}
 if only:
  results=json.loads((ROOT/'PHOTODISINTEGRATION_RESULTS.json').read_text())
  results['cases']=[x for x in results['cases'] if x['name']!=only]
 for name,tau,y,bree,braa,ne,nt in cases:
  if only and name!=only:continue
  # These modules import the settings by value. Set their effective settings,
  # rather than only changing params after the imports have already occurred.
  params.NE_pd=cascade.NE_pd=ne;params.NT_pd=nucl.NT_pd=nt
  assert cascade.NE_pd==ne and nucl.NT_pd==nt
  n0a=7.04*y # post-electron-annihilation entropy/photon number at T0=10 keV
  print('ACROPOLIS starting',name,'tau',tau,'n_phi/n_gamma',n0a,flush=True)
  m=DecayModel(1000*mass,tau,.01,n0a,bree,braa)
  covered=(m._sII.temperature_range()[0]<=m._sTrg[0]<=m._sTrg[1]<=m._sII.temperature_range()[1])
  if not covered:raise ValueError('Decay outside background temperature table')
  t0=time.monotonic();final=m.run_disintegration();duration=time.monotonic()-t0
  baseline=m._squeeze_decays(m._sII.bbn_abundances())
  bh=baseline[:,0];fh=final[:,0]
  change={'D_over_H_relative':float((fh[2]/fh[1])/(bh[2]/bh[1])-1),
          'He4_mass_fraction_relative':float(fh[5]/bh[5]-1),
          'He3_over_D_relative':float((fh[4]/fh[2])/(bh[4]/bh[2])-1)}
  row={'name':name,'lifetime_seconds':tau,'yield_per_entropy':y,'n0_over_photons_at_10keV':n0a,
       'branch_electrons':bree,'branch_photons':braa,'energy_grid_points_per_decade':ne,
       'temperature_grid_points_per_decade':nt,'effective_grid_settings_verified':True,
       'covered_temperature_range_MeV':list(m._sTrg),
       'baseline':baseline.tolist(),'final':final.tolist(),'relative_changes':change,
       'runtime_seconds_non_scientific':duration}
  results['cases'].append(row)
  (ROOT/'PHOTODISINTEGRATION_RESULTS.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
  print('ACROPOLIS finished',name,change,'seconds',duration,flush=True)
 for module in ['models.py','nucl.py','cascade.py','params.py','obs.py']:
  f=__import__('pathlib').Path(__import__('acropolis').__file__).parent/module
  if f.exists():results.setdefault('source_sha256',{})[module]=hashlib.sha256(f.read_bytes()).hexdigest()
 (ROOT/'PHOTODISINTEGRATION_RESULTS.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--only');run(p.parse_args().only)
