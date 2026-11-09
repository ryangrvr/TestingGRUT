"""Independent archived-data audit: no production imports or source evaluation."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import numpy as np

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root',type=Path,required=True)
parser.add_argument('--raw',type=Path,required=True)
parser.add_argument('--r1-point',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
FROZEN=args.source_root
ROOT=FROZEN/'research/reality_framework/grut_working_theory_20261109'
OUT=args.raw
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def complex_array(encoded):
    a=np.asarray(encoded,dtype=float)
    return a[...,0]+1j*a[...,1]
def exact(record):return F(record['numerator'],record['denominator'])
def rational_atan_bounds(x):
    partial=F(0)
    for k in range(16):partial+=((-1)**k)*x**(2*k+1)/F(2*k+1)
    next_partial=partial+x**33/F(33)
    return partial,next_partial
def maximum_difference(a,b):return float(np.max(abs(np.asarray(a)-np.asarray(b))))

d=read(OUT/'spectral.json');cfg=read(ROOT/'spectral_protocol.json')
provenance=read(OUT/'provenance.json');resource=read(OUT/'spectral_resources.json')
r1=read(args.r1_point);r1ref=r1['refined_reference']
slack=cfg['gates']['floating_comparison_slack'];consistency=[]
source_hashes={path:hashlib.sha256((FROZEN/path).read_bytes()).hexdigest()==value
 for path,value in provenance['source_sha256'].items()
 if 'grut_working_theory_20261109/spectral_' in path or path.endswith('/r1_reference.py')}

# Independent exact arithmetic; a narrower sixteen-term Machin bracket
# verifies inclusion in the archived eight-term bracket.
lo5,hi5=rational_atan_bounds(F(1,5));lo239,hi239=rational_atan_bounds(F(1,239))
pi_lo=16*lo5-4*hi239;pi_hi=16*hi5-4*lo239
cert=d['rational_certificate'];mstar=F(1323,163328);qstar=F(1,16)
r2=14*(F(1,23**2)+4*F(22,7)**2/F(3*23**3))
fraction_checks={
 'strict_pi_enclosure':3<exact(cert['pi_lower'])<=pi_lo<pi_hi<=exact(cert['pi_upper'])<F(22,7),
 'coefficient_loss':exact(cert['coefficient_loss_squared_upper'])==r2 and r2<F(1,16),
 'covariance_lower':exact(cert['covariance_lower_bound'])==F(9,40)*F(9,16)/(2*F(22,7)*F(87,35))==mstar,
 'covariance_upper':exact(cert['covariance_upper_bound'])==qstar}
for count in (16,64,256,1024):
    row=cert['mode_bounds'][str(count)];tail=F(7,9*count*count)
    bound=tail*(1/mstar+qstar/mstar**2)
    fraction_checks[f'mode_{count}_exact_tail']=exact(row['raw_primary_tail'])==tail
    fraction_checks[f'mode_{count}_exact_bound']=exact(row['uniform_canonical_error'])==bound
    fraction_checks[f'mode_{count}_classification']=row['below_one_tenth_percent']==(bound<F(1,1000))
    consistency.extend([abs(float(tail)-row['raw_primary_tail']['decimal']),abs(float(bound)-row['uniform_canonical_error']['decimal'])])
primary_bound=F(7,9*1024**2)*(1/mstar+qstar/mstar**2)
fraction_checks['primary_exact_claim']=exact(d['primary_mathematical_uniform_certificate']['bound'])==primary_bound<F(1,1000)

cases={**d['predictors'],'reference':d['reference'],'refined_reference':d['refined_reference']}
computed={};case_checks={}
for name,case in cases.items():
    C=np.asarray(case['covariance'],dtype=float);values,U=np.linalg.eigh(C[:3,:3])
    W=U@np.diag(values**-.5)@U.T
    R={t:W@complex_array(M)[:3,:3]@W for t,M in case['responses'].items()}
    scalar={t:np.diag(complex_array(M))/np.diag(C) for t,M in case['responses'].items()}
    weights=np.asarray(case['source_weights']);edges=np.asarray(case['source_boundaries']);a=5/23
    height=weights/a
    TV=abs(height[:,0])+abs(height[:,-1])+np.sum(abs(np.diff(height,axis=1)),axis=1)
    T=TV*TV*25/(2*np.pi**4*case['modes']**2);E=float(T[:3].sum())
    canonical=E/values[0]+E*values[-1]/values[0]**2
    entry=np.sqrt(T[:,None]*T[None,:]);scalar_tail=2*T/np.diag(C)
    t0=float(np.linalg.norm(R['0']-np.eye(3),2))
    contraction=max(0.,max(float(np.linalg.norm(M,2))-1 for M in R.values()))
    weights_error=maximum_difference(weights,r1['probe_weights'])
    boundaries_error=maximum_difference(edges,r1ref['probe_boundary_integrals'])
    constants_error=maximum_difference(case['source_constants'],r1ref['probe_constants'])
    partition=abs(weights.sum(axis=1)+edges.sum(axis=1)-[1.,0.,1/3,1.])
    saved_errors=[maximum_difference(values,case['primary']['eigenvalues']),
      abs(values[-1]/values[0]-case['primary']['condition']),
      maximum_difference(T,case['variance_tail_bounds']),maximum_difference(entry,case['entry_tail_bounds']),
      abs(E-case['raw_primary_tail_bound']),abs(canonical-case['uniform_canonical_tail_bound']),
      maximum_difference(scalar_tail,case['uniform_scalar_tail_bounds']),maximum_difference(TV,case['source_total_variations']),
      maximum_difference(partition,case['source_partition_errors']),abs(t0-case['primary']['t0_identity_error']),
      abs(contraction-case['primary']['contraction_excess'])]
    for t in R:
        saved_errors.extend([maximum_difference(R[t],complex_array(case['primary']['responses'][t])),
          maximum_difference(scalar[t],complex_array(case['scalar_responses'][t]))])
    consistency.extend(saved_errors)
    checks={'positive_primary_Gram':bool(values[0]>cfg['gram_minimum_eigenvalue']),
      'condition_limit':bool(values[-1]/values[0]<cfg['gram_condition_limit']),
      'raw_t0_equals_covariance':maximum_difference(complex_array(case['responses']['0']),C)<=slack,
      'covariance_symmetric':maximum_difference(C,C.T)<=slack,
      'canonical_t0_identity':t0<=cfg['gates']['t0_identity'],
      'contraction':contraction<=cfg['gates']['contraction'],
      'apriori_Gram_envelopes':bool(values[0]>=float(mstar)-slack and values[-1]<=float(qstar)+slack),
      'apriori_variation_envelopes':bool(np.all(TV[:3]<=np.arange(1,4)*np.pi/5+slack)),
      'source_partition':bool(max(partition)<=slack),
      'archived_R1_source_dictionary_unchanged':max(weights_error,boundaries_error,constants_error)<=slack,
      'saved_whitening_and_tail_metrics_reproduce':max(saved_errors)<=slack}
    case_checks[name]={'checks':checks,'saved_metric_max_error':max(saved_errors),
      'R1_weight_error':weights_error,'R1_boundary_error':boundaries_error,'R1_constant_error':constants_error,
      'recomputed_canonical_tail_bound':canonical,'recomputed_raw_primary_tail_bound':E,'passed':all(checks.values())}
    computed[name]={'C':C,'R':R,'scalar':scalar,'T':T,'E':E,'canonical':canonical,'entry':entry,'scalar_tail':scalar_tail}

def audit_comparison(label,left_name,right_name,stored):
    left=cases[left_name];right=cases[right_name];A=computed[left_name];B=computed[right_name]
    raw_allow=A['E']+B['E'];canon_allow=A['canonical']+B['canonical']
    entry_allow=A['entry']+B['entry'];scalar_allow=A['scalar_tail']+B['scalar_tail']
    checks={};rawmax=canonmax=scalarmax=0.;local_errors=[]
    for t in left['responses']:
        delta=complex_array(left['responses'][t])-complex_array(right['responses'][t])
        white=A['R'][t]-B['R'][t];scalar=A['scalar'][t]-B['scalar'][t]
        raw=float(np.linalg.norm(delta[:3,:3],2));canon=float(np.linalg.norm(white,2));se=abs(scalar)
        rawmax=max(rawmax,raw);canonmax=max(canonmax,canon);scalarmax=max(scalarmax,float(max(se)))
        accepted=bool(raw<=raw_allow+slack and canon<=canon_allow+slack
          and np.all(abs(delta)<=entry_allow+slack) and np.all(se<=scalar_allow+slack))
        row=stored['times'][t]
        local_errors.extend([abs(raw-row['raw_primary_operator_error']),abs(canon-row['canonical_primary_operator_error']),
          maximum_difference(abs(delta),row['entry_absolute_errors']),maximum_difference(se,row['scalar_absolute_errors']),
          maximum_difference(white,complex_array(row['canonical_signed_difference'])),
          maximum_difference(scalar,complex_array(row['scalar_signed_difference']))])
        checks[t]={'raw_operator_error':raw,'canonical_operator_error':canon,
          'max_scalar_error':float(max(se)),'bound_compliance':accepted,
          'stored_classification_matches':accepted==row['within_tail_allowances']}
    local_errors.extend([abs(rawmax-stored['max_raw_primary_error']),abs(canonmax-stored['max_canonical_primary_error']),
      abs(scalarmax-stored['max_scalar_error']),abs(raw_allow-stored['raw_primary_allowance']),
      abs(canon_allow-stored['canonical_primary_allowance']),maximum_difference(entry_allow,stored['entry_allowances']),
      maximum_difference(scalar_allow,stored['scalar_allowances'])])
    consistency.extend(local_errors)
    all_bounds=all(v['bound_compliance'] for v in checks.values())
    return {'times':checks,'maximum_saved_metric_error':max(local_errors),'maximum_canonical_error':canonmax,
      'maximum_scalar_error':scalarmax,'canonical_predictor_plus_reference_allowance':canon_allow,
      'bound_classification_matches':all_bounds==stored['bound_compliance_passed'],
      'passed':bool(all_bounds and all(v['stored_classification_matches'] for v in checks.values()) and max(local_errors)<=slack)}

comparisons={name:audit_comparison(name,name,'refined_reference',d['comparisons'][name]) for name in d['predictors']}
doubling=audit_comparison('reference_doubling','reference','refined_reference',d['reference_doubling'])
pre=read(ROOT/'spectral_precheck.json')
quality={'synthetic_checks':bool(pre['precheck']['passed'] and all(pre['precheck']['checks'].values())),
 'exact_rational_inequalities':all(fraction_checks.values()),
 'source_partition':all(r['checks']['source_partition'] for r in case_checks.values()),
 't0_normalization':all(r['checks']['canonical_t0_identity'] for r in case_checks.values()),
 'response_contraction':all(r['checks']['contraction'] for r in case_checks.values()),
 'analytic_Gram_envelopes':all(r['checks']['apriori_Gram_envelopes'] for r in case_checks.values()),
 'analytic_source_variation_envelopes':all(r['checks']['apriori_variation_envelopes'] for r in case_checks.values()),
 'reference_doubling':max(doubling['maximum_canonical_error'],doubling['maximum_scalar_error'])<=cfg['gates']['reference_doubling'],
 'independent_R1_implementation':d['independent_R1_check']['maximum']<=cfg['gates']['independent_R1_agreement'],
 'tail_bound_compliance':all(r['passed'] for r in comparisons.values()) and doubling['passed']}
quality_matches=quality==d['quality_gates']
registered=cfg['N']==22 and cfg['length']==5 and cfg['g']==1 and cfg['mass_R']==0 and cfg['all_predictor_modes']==[16,64,256,1024]
registered=registered and d['primary_modes_fixed_before_evaluation']==1024 and cfg['response_times']==[0,.17,.61,1.09,2.03,3.17]
registered=registered and cases['reference']['modes']==16384 and cases['refined_reference']['modes']==32768
scope=all(d[k] is False for k in ('full_physical_dictionary_proved','lattice_continuum_mismatch_repaired',
 'nonlinear_detector_instrument_validated','empirical_theory_validated'))
result={'scope':'Independent postrun audit of archived matrices and exact rational arithmetic only. No production module imported; no source coefficients, continuum evolution, or lattice evolution rerun.',
 'freeze_commit':provenance['freeze_commit'],'spectral_sha256':hashlib.sha256((OUT/'spectral.json').read_bytes()).hexdigest(),
 'registered_configuration':bool(registered),'source_hash_checks':source_hashes,'exact_fraction_checks':fraction_checks,
 'case_audits':case_checks,'comparison_audits':comparisons,'reference_doubling_audit':doubling,
 'reconstructed_quality_gates':quality,'all_ten_quality_classifications_match':quality_matches,
 'maximum_recomputed_metric_difference':max(consistency),
 'primary_measured_error_at_six_times':comparisons['1024']['maximum_canonical_error'],
 'primary_exact_all_time_truncation_bound':{'numerator':primary_bound.numerator,'denominator':primary_bound.denominator,'decimal':float(primary_bound)},
 'resource_record_within_limits':resource['returncode']==0 and resource['wall_seconds']<=120 and resource['peak_private_bytes']<=1073741824,
 'limitations_preserved':scope,
 'time_and_precision_limit':'The analytic modal-truncation theorem covers all real times for exact mathematical phases and canonical roots. Executable results were checked only at six fixed times0,.17,.61,1.09,2.03,3.17. Their floating error is not interval-certified uniformly in time; phase multiplication/argument reduction can accumulate finite-precision error at very large times. No uniform all-time machine-arithmetic certificate is claimed.',
 'independence_limit':'The frozen R1 raw-response agreement and synthetic precheck are archived cross-checks inspected here; this audit independently reconstructs all W1 whitening, comparison, and tail metrics but does not rerun either source implementation.',
 'passed':bool(registered and all(source_hashes.values()) and all(fraction_checks.values())
   and all(r['passed'] for r in case_checks.values()) and all(r['passed'] for r in comparisons.values())
   and doubling['passed'] and quality_matches and all(quality.values()) and scope and max(consistency)<=slack)}
path=args.output
with path.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps({k:result[k] for k in ('passed','registered_configuration','all_ten_quality_classifications_match',
 'maximum_recomputed_metric_difference','primary_measured_error_at_six_times','primary_exact_all_time_truncation_bound','limitations_preserved')},indent=2))
