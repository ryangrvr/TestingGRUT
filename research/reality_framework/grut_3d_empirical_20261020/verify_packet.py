"""Integrity/reproduction checks, distinct from scientific pass/fail criteria."""
from pathlib import Path
import hashlib, json
from unittest.mock import patch
import numpy as np
import empirical, empirical_gravity
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    spec=json.loads((P/'protocol.json').read_text());trans=json.loads((P/'TRANSVERSE_PROTOCOL.json').read_text())
    checks={name:sha(P/name)==v for name,v in spec['files_sha256'].items()}
    for name,prefix in [('SPECTRAL_PROTOCOL.json','spectral.py'),('LARGE_DOMAIN_PROTOCOL.json','large_domain_modes.cjs'),('TRANSVERSE_PROTOCOL.json','run_transverse.py')]:checks[prefix]=sha(P/prefix)==json.loads((P/name).read_text())['source_sha256']
    cases=spec['cases']+trans['cases']
    for c in cases:
      p=P/'raw'/(c['id']+'.csv');a=np.genfromtxt(p,delimiter=',',names=True,ndmin=1);rec=json.loads((P/'raw'/(c['id']+'.execution.json')).read_text())
      checks[c['id']+'_completed']=bool(rec['returncode']==0 and abs(a[-1]['t']-c['T'])<1e-10 and np.isfinite(np.array(a.tolist(),dtype=float)).all())
      checks[c['id']+'_source']=rec['source_sha256']==spec['files_sha256']['cartesian.cpp']
    original=json.loads((P/'LARGE_DOMAIN_PROTOCOL.json').read_text())
    parent=P/'../grut_candidate1_two_object_20261019/persistence/kernel_modes.cjs'
    checks['original_kernel_preserved']=sha(parent)==original['parent_sha256']
    for h in original['grids']:
      d=json.loads((P/'raw'/('large_modes_h'+format(h,'.6f')+'.json')).read_text())
      checks['large_domain_'+str(h)]=d['trace'][-1]['time']==original['T'] and d['parent_sha256']==original['parent_sha256']
    # Restarted fine trajectory agrees with the preserved pre-disconnection prefix.
    partial=np.genfromtxt(P/'raw/interrupted_pair_reference_fine.csv',delimiter=',',names=True)
    full=np.genfromtxt(P/'raw/pair_reference_fine.csv',delimiter=',',names=True);ix=np.searchsorted(full['t'],partial['t'])
    physical=['E','Q','c0_Q12','c1_Q12','tan0_norm','tan1_norm','tan2_norm']
    error=max(float(np.max(abs(partial[k]-full[k][ix]))/max(np.max(abs(full[k][ix])),1e-30)) for k in physical)
    checks['interruption_replay']=error<1e-11
    s=json.loads((P/'SPECTRAL_RESULTS.json').read_text());v=s['radial_spectra']
    ratios={name:[abs(v[0][name][0]/v[1][name][0]),abs(v[1][name][0]/v[2][name][0])] for name in ['phase_l0','translation_l1']}
    checks['phase_zero_second_order']=all(3.9<x<4.1 for x in ratios['phase_l0'])
    checks['translation_zero_second_order']=all(3.9<x<4.1 for x in ratios['translation_l1'])
    checks['spectral_fixed_charge_positive']=all(r['fixed_charge_l0'][0]>0 for r in v)
    checks['spectral_shape_positive']=all(r['shape_l2'][0]>0 for r in v)
    maxeig=max(max(r[k]) for r in v for k in r if k.endswith('_residuals'))
    checks['eigenpair_residuals']=maxeig<1e-7
    # Enforce operational separation from held-out data in the prediction code.
    real_read=Path.read_text
    def guarded(self,*a,**kw):
      if self.name.startswith('heldout_'):raise RuntimeError('prediction read held-out data')
      return real_read(self,*a,**kw)
    with patch.object(Path,'read_text',guarded):
      ep=empirical.predict();gp=empirical_gravity.prediction()
    checks['prediction_does_not_read_heldout']=bool(ep['omega_hidden_floor']>0 and gp['H0_lower_bound_km_s_Mpc']>0)
    if (P/'RAW_MANIFEST.json').exists():
      raw=json.loads((P/'RAW_MANIFEST.json').read_text());checks['raw_manifest']=all(sha(P/'raw'/n)==r['sha256'] for n,r in raw['raw_files'].items())
    out={'checks':checks,'all_integrity_checks_pass':all(checks.values()),'count':len(checks),'interruption_maximum_relative_prefix_difference':error,'zero_mode_reduction_factors':ratios,'maximum_eigenpair_absolute_residual':maxeig,'scope':'These are integrity, completed-execution, operator and prediction-input checks. A passing check does not turn failed symmetry-convergence science criteria into successes.'}
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
    if not all(checks.values()):raise RuntimeError('verification failed')
if __name__=='__main__':main()
