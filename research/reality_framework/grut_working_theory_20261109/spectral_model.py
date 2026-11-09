"""W1 frozen-family spectral predictor; physical evaluation is explicit only.

Running this file performs synthetic checks only. evaluate_study() is the
separate, post-freeze physical-source entry point. No lattice states load.
"""
import argparse, importlib.util, json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

NAMES=('sine1','sine2','sine3','beta_left')

def encode(x):
    x=np.asarray(x,dtype=complex)
    return np.stack((x.real,x.imag),axis=-1).tolist()

def decode(x):
    a=np.asarray(x)
    return a[...,0]+1j*a[...,1] if a.ndim==3 and a.shape[-1]==2 else np.asarray(x,dtype=complex)

def fraction_record(x):
    return {'numerator':x.numerator,'denominator':x.denominator,'decimal':float(x)}

def atan_bracket(x,terms=8):
    """Exact alternating-series bracket; x is a rational in (0,1)."""
    total=sum(((-1)**k)*x**(2*k+1)/Q(2*k+1) for k in range(terms))
    following=total+((-1)**terms)*x**(2*terms+1)/Q(2*terms+1)
    return min(total,following),max(total,following)

def rational_certificate(modes=(16,64,256,1024)):
    """Analytic N22/L5 envelope only; never evaluates source coefficients."""
    l5,u5=atan_bracket(Q(1,5));l239,u239=atan_bracket(Q(1,239))
    lower_pi=16*l5-4*u239;upper_pi=16*u5-4*l239
    r2=14*(Q(1,23**2)+Q(4,3)*Q(22,7)**2/Q(23**3))
    minimum=Q(9,40)*Q(9,16)/(2*Q(22,7)*Q(87,35));maximum=Q(1,16)
    bounds={}
    for m in modes:
        E=Q(7,9*m*m);B=E*(1/minimum+maximum/minimum**2)
        bounds[str(m)]={'raw_primary_tail':fraction_record(E),
          'uniform_canonical_error':fraction_record(B),
          'below_one_tenth_percent':bool(B<Q(1,1000))}
    proof={'machin_lower_above_three':lower_pi>3,'machin_upper_below_22_over_7':upper_pi<Q(22,7),
      'coefficient_loss_squared_below_one_sixteenth':r2<Q(1,16),
      'positive_covariance_lower_bound':minimum>0,
      'primary_1024_below_one_thousandth':bounds['1024']['below_one_tenth_percent']}
    return {'scope':'Exact rational verification of analytic envelopes; no source response evaluated.',
      'pi_lower':fraction_record(lower_pi),'pi_upper':fraction_record(upper_pi),
      'coefficient_loss_squared_upper':fraction_record(r2),
      'covariance_lower_bound':fraction_record(minimum),'covariance_upper_bound':fraction_record(maximum),
      'mode_bounds':bounds,'proof_checks':proof,'passed':all(proof.values())}

def native_sources(length,n):
    """Independent closed-form cell integrals; affine edge terms preserved."""
    a=length/(n+1);edges=np.arange(1,n+1,dtype=float)*a
    W=[];edge=[]
    for k in (1,2,3):
        ang=k*np.pi/length
        W.append((np.cos(ang*edges[:-1])-np.cos(ang*edges[1:]))/(2*k))
        edge.append([(1-np.cos(ang*a))/(2*k),(np.cos(ang*(length-a))-(-1)**k)/(2*k)])
    def primitive(z):
        x=np.asarray(z)/length
        return 15*x*x-40*x**3+45*x**4-24*x**5+5*x**6
    W.append(np.diff(primitive(edges)));edge.append([primitive(a),1-primitive(length-a)])
    W=np.asarray(W);edge=np.asarray(edge);height=W/a
    variation=abs(height[:,0])+np.sum(abs(np.diff(height,axis=1)),axis=1)+abs(height[:,-1])
    integral=np.array([1.,0.,1/3,1.])
    return {'a':a,'edges':edges,'weights':W,'height':height,'boundaries':edge,
      'constants':edge.sum(axis=1)/4,'variation':variation,
      'partition_errors':abs(W.sum(axis=1)+edge.sum(axis=1)-integral)}

def whitened_response(covariance,responses,floor=1e-10,condition_limit=1e6):
    C=np.asarray(covariance,dtype=float);e,U=np.linalg.eigh((C+C.T)/2)
    if e[0]<=floor or e[-1]/e[0]>=condition_limit:raise ValueError('Frozen primary Gram condition failed')
    inverse=(U*(1/np.sqrt(e)))@U.T
    R={t:inverse@decode(value)@inverse for t,value in responses.items()}
    return {'eigenvalues':e.tolist(),'condition':float(e[-1]/e[0]),
      'responses':{t:encode(value) for t,value in R.items()},
      'contraction_excess':max(0.,max(float(np.linalg.norm(value,2))-1 for value in R.values())),
      't0_identity_error':float(np.linalg.norm(R['0']-np.eye(C.shape[0]),2))}

def compute_case(card,modes):
    """Physical-source evaluation; call only after the protocol freeze."""
    if card['N']!=22 or card['length']!=5 or card['g']!=1 or card['mass_R']!=0:
        raise ValueError('Unregistered source model')
    length=float(card['length']);n=card['N'];data=native_sources(length,n)
    wave=np.arange(1,modes+1,dtype=float)*np.pi/length
    omega=np.sqrt(wave*wave+1/np.pi)
    edges=data['edges'];cell_integrals=np.sqrt(2/length)*(np.cos(edges[:-1,None]*wave)-np.cos(edges[1:,None]*wave))/wave
    coefficients=data['height']@cell_integrals
    weighted=coefficients/np.sqrt(2*np.pi*omega)
    covariance=weighted@weighted.T
    responses={f'{float(t):g}':encode((weighted*np.exp(-1j*omega*t))@weighted.T) for t in card['response_times']}
    primary={t:encode(decode(value)[:3,:3]) for t,value in responses.items()}
    white=whitened_response(covariance[:3,:3],primary,card['gram_minimum_eigenvalue'],card['gram_condition_limit'])
    tails=data['variation']**2*length**2/(2*np.pi**4*modes*modes)
    entry=np.sqrt(tails[:,None]*tails[None,:]);E=float(tails[:3].sum())
    m=white['eigenvalues'][0];q=white['eigenvalues'][-1]
    variance=np.diag(covariance)
    if np.any(variance<=0):raise ValueError('Nonpositive retained scalar variance')
    scalar={t:encode(np.diag(decode(value))/variance) for t,value in responses.items()}
    return {'modes':modes,'source_weights':data['weights'].tolist(),'source_boundaries':data['boundaries'].tolist(),
      'source_constants':data['constants'].tolist(),'source_partition_errors':data['partition_errors'].tolist(),
      'source_total_variations':data['variation'].tolist(),'leading_frequencies':omega[:8].tolist(),
      'leading_coefficients':coefficients[:,:8].tolist(),'covariance':covariance.tolist(),'responses':responses,
      'primary':white,'scalar_responses':scalar,'variance_tail_bounds':tails.tolist(),
      'entry_tail_bounds':entry.tolist(),'raw_primary_tail_bound':E,
      'uniform_canonical_tail_bound':float(E/m+E*q/(m*m)),
      'uniform_scalar_tail_bounds':(2*tails/variance).tolist(),
      'floating_scope':'Values and data-dependent tail estimates use ordinary floating arithmetic; the separate rational theorem certifies the exact mathematical family.'}

def decode_vector(value):
    a=np.asarray(value);return a[:,0]+1j*a[:,1]

def compare_cases(left,right,slack):
    raw_bound=left['raw_primary_tail_bound']+right['raw_primary_tail_bound']
    canonical_bound=left['uniform_canonical_tail_bound']+right['uniform_canonical_tail_bound']
    entry_bound=np.asarray(left['entry_tail_bounds'])+right['entry_tail_bounds']
    scalar_bound=np.asarray(left['uniform_scalar_tail_bounds'])+right['uniform_scalar_tail_bounds']
    records={};passed=True
    for t in left['responses']:
        raw=decode(left['responses'][t])-decode(right['responses'][t])
        canon=decode(left['primary']['responses'][t])-decode(right['primary']['responses'][t])
        scalar=decode_vector(left['scalar_responses'][t])-decode_vector(right['scalar_responses'][t])
        raw_error=float(np.linalg.norm(raw[:3,:3],2));canon_error=float(np.linalg.norm(canon,2))
        scalar_error=abs(scalar);entry_error=abs(raw)
        ok=bool(raw_error<=raw_bound+slack and canon_error<=canonical_bound+slack
          and np.all(entry_error<=entry_bound+slack) and np.all(scalar_error<=scalar_bound+slack))
        passed=passed and ok
        records[t]={'raw_primary_operator_error':raw_error,'canonical_primary_operator_error':canon_error,
          'entry_absolute_errors':entry_error.tolist(),'scalar_absolute_errors':scalar_error.tolist(),
          'canonical_signed_difference':encode(canon),'scalar_signed_difference':encode(scalar),'within_tail_allowances':ok}
    return {'times':records,'raw_primary_allowance':raw_bound,'canonical_primary_allowance':canonical_bound,
      'entry_allowances':entry_bound.tolist(),'scalar_allowances':scalar_bound.tolist(),'slack':slack,
      'max_raw_primary_error':max(x['raw_primary_operator_error'] for x in records.values()),
      'max_canonical_primary_error':max(x['canonical_primary_operator_error'] for x in records.values()),
      'max_scalar_error':max(v for x in records.values() for v in x['scalar_absolute_errors']),
      'bound_compliance_passed':passed}

def independent_R1_check(card,case,path):
    spec=importlib.util.spec_from_file_location('W1_frozen_R1_reference',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    other=module.reference(card,card['N'],case['modes'])
    differences={'weights':float(np.max(abs(np.asarray(case['source_weights'])-other['probe_weights']))),
      'covariance':float(np.max(abs(np.asarray(case['covariance'])-other['covariance']))),
      'raw_responses':max(float(np.max(abs(decode(case['responses'][t])-decode(other['responses'][t])))) for t in case['responses']),
      'canonical_responses':max(float(np.max(abs(decode(case['primary']['responses'][t])-decode(other['primary_diagnostics']['whitened_responses'][t])))) for t in case['responses'])}
    return {'reference_file':str(path),'errors':differences,'maximum':max(differences.values()),
      'scope':'Independent frozen R1 implementation sanity check; it is not the uniform-error proof.'}

def evaluate_study(card,r1_reference_path):
    """The sole post-freeze bounded physical entry point; no writes or reruns."""
    expected=[16,64,256,1024]
    if card['all_predictor_modes']!=expected or card['primary_modes']!=1024:
        raise ValueError('Primary mode count must remain fixed prospectively')
    if card['response_times']!=[0,.17,.61,1.09,2.03,3.17] or card['reference_modes']!=16384 or card['refined_reference_modes']!=32768:
        raise ValueError('Unregistered W1 times or reference counts')
    certificate=rational_certificate(expected);pre=synthetic_precheck()
    cases={str(m):compute_case(card,m) for m in expected}
    reference=compute_case(card,card['reference_modes']);fine=compute_case(card,card['refined_reference_modes'])
    slack=card['gates']['floating_comparison_slack']
    comparisons={key:compare_cases(value,fine,slack) for key,value in cases.items()}
    doubling=compare_cases(reference,fine,slack)
    independent=independent_R1_check(card,reference,r1_reference_path)
    allcases=list(cases.values())+[reference,fine]
    mstar=certificate['covariance_lower_bound']['decimal'];qstar=certificate['covariance_upper_bound']['decimal']
    quality={'synthetic_checks':pre['passed'],'exact_rational_inequalities':certificate['passed'],
      'source_partition':all(max(c['source_partition_errors'])<=slack for c in allcases),
      't0_normalization':all(c['primary']['t0_identity_error']<=card['gates']['t0_identity'] for c in allcases),
      'response_contraction':all(c['primary']['contraction_excess']<=card['gates']['contraction'] for c in allcases),
      'analytic_Gram_envelopes':all(c['primary']['eigenvalues'][0]>=mstar-slack and c['primary']['eigenvalues'][-1]<=qstar+slack for c in allcases),
      'analytic_source_variation_envelopes':all(np.all(np.asarray(c['source_total_variations'][:3])<=np.arange(1,4)*np.pi/card['length']+slack) for c in allcases),
      'reference_doubling':max(doubling['max_canonical_primary_error'],doubling['max_scalar_error'])<=card['gates']['reference_doubling'],
      'independent_R1_implementation':independent['maximum']<=card['gates']['independent_R1_agreement'],
      'tail_bound_compliance':all(c['bound_compliance_passed'] for c in comparisons.values()) and doubling['bound_compliance_passed']}
    bound=certificate['mode_bounds']['1024']['uniform_canonical_error']
    exact_claim=Q(bound['numerator'],bound['denominator'])<Q(1,1000)
    return {'id':card['id'],'complete':True,'numerically_valid':all(quality.values()),'quality_gates':quality,
      'primary_modes_fixed_before_evaluation':1024,'primary_mathematical_uniform_certificate':{
        'passed':bool(exact_claim),'threshold':.001,'all_real_times':True,'bound':bound,
        'scope':'Exact-arithmetic mathematical modal truncation theorem, not interval certification of executable floating arithmetic or physical dictionary validity.'},
      'rational_certificate':certificate,'predictors':cases,'reference':reference,'refined_reference':fine,
      'comparisons':comparisons,'reference_doubling':doubling,'independent_R1_check':independent,
      'verdict':'Working finite-mode predictor of the specified continuum source response only.' if exact_claim and all(quality.values()) else 'Report mathematical and implementation outcomes separately; do not tune or extend.',
      'full_physical_dictionary_proved':False,'lattice_continuum_mismatch_repaired':False,
      'nonlinear_detector_instrument_validated':False,'empirical_theory_validated':False}

def synthetic_precheck():
    from scipy.linalg import expm
    frequency=np.array([.7,1.3,2.1,4.2]);A=np.array([[.3,.08,-.03,.02],[.05,.27,.04,-.03],[-.02,.06,.22,.01]])
    rotation=np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]],dtype=float)/2
    H=rotation@np.diag(frequency)@rotation.T;vectors=rotation@A.T
    times=[-1.17,0.,.43,2.71]
    full={f'{t:g}':encode((A*np.exp(-1j*frequency*t))@A.T) for t in times}
    retained=A[:,:3];small={f'{t:g}':encode((retained*np.exp(-1j*frequency[:3]*t))@retained.T) for t in times}
    C=A@A.T;CM=retained@retained.T;Delta=C-CM;E=float(np.linalg.norm(Delta,2))
    wf=whitened_response(C,full);wm=whitened_response(CM,small)
    m,q=wm['eigenvalues'][0],wm['eigenvalues'][-1];bound=E/m+E*q/(m*m)
    dense_error=max(float(np.linalg.norm(vectors.T@expm(-1j*H*t)@vectors-decode(full[f'{t:g}']),2)) for t in times)
    raw=max(float(np.linalg.norm(decode(full[k])-decode(small[k]),2)) for k in full)
    white=max(float(np.linalg.norm(decode(wf['responses'][k])-decode(wm['responses'][k]),2)) for k in full)
    checks={'dense_independent_evolution':dense_error<1e-12,'discarded_covariance_PSD':np.linalg.eigvalsh(Delta)[0]>-1e-14,
      'raw_uniform_tail_example':raw<=E+1e-14,'canonical_tail_example':white<=bound+1e-14,
      't0_preserved':max(wf['t0_identity_error'],wm['t0_identity_error'])<1e-12,
      'contraction':max(wf['contraction_excess'],wm['contraction_excess'])<1e-12,
      'rational_primary_certificate':rational_certificate()['passed']}
    checks={key:bool(value) for key,value in checks.items()}
    return {'scope':'Synthetic four-mode matrices and exact rational analytic envelopes only; no physical source evaluated.',
      'checks':checks,'dense_error':dense_error,'raw_error':raw,'raw_bound':E,
      'canonical_error':white,'canonical_bound':bound,'passed':all(checks.values())}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--synthetic-output',type=Path,required=True)
    args=parser.parse_args();result={'precheck':synthetic_precheck(),'rational_certificate':rational_certificate()}
    with args.synthetic_output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result,indent=2))
