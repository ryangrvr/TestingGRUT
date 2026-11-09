"""Postrun algebra only, using the archived registered N22 current matrices.

No new source strength is simulated; no frozen file is changed.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--raw',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
PRODUCTION=args.raw
def relative_norm(matrix,reference):
    w,U=np.linalg.eigh((reference+reference.T)/2)
    if w[0]<=0:raise ValueError('Nonpositive primary current matrix')
    inverse=(U*(1/np.sqrt(w)))@U.T
    return float(np.linalg.norm(inverse@(matrix-reference)@inverse,2))

def derive(label):
    path=PRODUCTION/(label+'.json');raw=path.read_bytes();d=json.loads(raw)
    weights=np.asarray(d['detector_weights'],float)
    hj=np.asarray(d['single_current_rules']['vacuum']['hopping_expectations'],float)
    v=np.asarray(d['probe_weights'],float)[:3]
    S0=np.asarray(d['single_current_rules']['vacuum']['fsum'],float)[:3,:3]
    Sc=np.asarray(d['refined_reference']['fsum'],float)[:3,:3]
    gram0=np.linalg.eigvalsh(S0);gramc=np.linalg.eigvalsh(Sc)
    e0=relative_norm(S0,Sc);wmax=float(np.max(abs(weights)));tau=.01;guard=.0001
    allowed=min(tau-guard,(tau-guard-e0)/(1+e0))
    conditions={'nonpositive_initial_bond_hoppings':bool(np.all(hj<=0)),
      'positive_vacuum_current_matrix':bool(gram0[0]>0),
      'positive_reference_current_matrix':bool(gramc[0]>0),
      'baseline_guarded_allowance':bool(e0+guard<tau),
      'registered_kick_cosine_check':bool(d['kick_cosine_identity_error']<=1e-9),
      'registered_double_kick_cosine_check':bool(d['twice_kick_cosine_identity_error']<=1e-9)}
    if not all(conditions.values()):raise ValueError('Sufficient-condition assumptions not supported by archive')
    quadratic=math.sqrt(2*allowed)/wmax if wmax else None
    cosine=math.acos(1-allowed)/wmax if wmax else None
    scan={}
    for key,row in d['instant_kick_scan'].items():
        kappa=float(key);S=np.asarray(row['fsum'],float)[:3,:3]
        delta=float(np.max(1-np.cos(kappa*weights)))
        scan[key]={'delta_exact_trigonometric_bound':delta,
          'delta_quadratic_bound':kappa*kappa*wmax*wmax/2,
          'disturbance_upper_bound':delta,'absolute_mismatch_upper_bound':e0+(1+e0)*delta,
          'archived_disturbance':relative_norm(S,S0),'archived_absolute_mismatch':relative_norm(S,Sc),
          'both_supported_by_sufficient_bound':bool(delta+guard<=tau and e0+(1+e0)*delta+guard<=tau)}
    return {'source_file':path.name,'source_sha256':hashlib.sha256(raw).hexdigest(),
      'N':d['N'],'refined':d['refined'],'conditions':conditions,
      'minimum_initial_bond_hopping':float(hj.min()),'maximum_initial_bond_hopping':float(hj.max()),
      'minimum_vacuum_current_eigenvalue':float(gram0[0]),'minimum_reference_current_eigenvalue':float(gramc[0]),
      'current_reconstruction_error':float(np.max(abs(-(v*hj)@v.T-S0))),
      'kick_cosine_identity_error':d['kick_cosine_identity_error'],
      'twice_kick_cosine_identity_error':d['twice_kick_cosine_identity_error'],
      'baseline_whitened_mismatch':e0,'maximum_absolute_detector_weight':wmax,
      'threshold':tau,'numerical_guard':guard,'maximum_allowed_delta':allowed,
      'quadratic_sufficient_radius':quadratic,'cosine_sufficient_radius':cosine,
      'registered_strength_bound_comparisons':scan}

points={label:derive(label) for label in ('N22','N22_refined')}
out={'study':'R1-01 postrun algebraic supplement','points':points,
  'conservative_radius_both_runs_quadratic':min(p['quadratic_sufficient_radius'] for p in points.values()),
  'conservative_radius_both_runs_cosine':min(p['cosine_sufficient_radius'] for p in points.values()),
  'scope':'A sufficient finite-matrix current-only condition for the instantaneous unitary source kick. Derived after registered production from archived matrices, with unchanged one-percent allowance and 1e-4 guard. No new physical trial, fitted parameter, cutoff change, or rerun. A guarded numerical implication of the stored matrices, not an interval-certified exact-ground-state or empirical theory result.',
  'excluded_claims':['Selective detector outcomes','Nonselective instrument channel unless separately specified',
    'Delayed states','General state-current correspondence','Correct Hamiltonian propagation',
    'Continuum theorem','Massive transfer','Experimental validation']}
path=args.output
with path.open('x',encoding='utf-8') as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(out,indent=2))
