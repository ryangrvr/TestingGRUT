"""Frozen R1 numerical validity and separate finite diagnostic hypotheses.

No aggregate physical/theory pass is defined. ``pass`` means numerical validity
only; every scientific claim is a separately reported, guarded hypothesis.
"""
import json, math
from pathlib import Path
import numpy as np
from shared import card, HERE
from r1_reference import decode_matrix, encode_matrix, scalar_diagnostics

GATE_DEFAULTS = {
    'numerical':1e-9, 'refinement':1e-8, 'reference_doubling':1e-7,
    'reference_tail':1e-4, 'single_probability_floor':.01,
    'double_probability_floor':.001, 'numerical_guard':1e-4,
    'static_relative':.01, 'normalized_response_absolute':.005,
    'centroid_relative':.005, 'cosine_eigen_absolute':.005,
    'current_relative':.01, 'galerkin_residual':.01,
    'instant_moment_absolute':.002, 'instant_moment_reserve':1e-4}

def sections(d):
    return {**{'single/'+t:v for t,v in d['single'].items()},'double':d['double']}

def moment_values(d):
    out={}
    for section,block in sections(d).items():
        for name,probe in block['probes'].items():
            for state,record in probe['states'].items():
                for field in ('mean','second_moment'):
                    out[f'{section}/{name}/{state}/{field}']=record[field]
            for field in ('first_matrix','second_matrix'):
                for i,row in enumerate(probe[field]):
                    for j,v in enumerate(row):out[f'{section}/{name}/{field}/{i}/{j}']=complex(*v)
    for name,row in d['probe_covariance'].items():
        for other,value in row.items():out[f'vacuum_covariance/{name}/{other}']=value
    return out

def moment_errors(d,r):
    p,q=moment_values(d),moment_values(r)
    if set(p)!=set(q):raise ValueError('Instantaneous observable schema differs')
    return {k:float(abs(p[k]-q[k])) for k in p}

def instant_tail_control(ref,cfg):
    """D1 Gaussian tail control, locally implemented to avoid import shadowing."""
    VF=ref['detector']['variance'];TF=ref['detector']['tail_variance_bound'];result={};maximum=0.
    for section,block in sections(ref).items():
        Q=(4 if section=='double' else 2)*cfg['coupling'];D=.5*Q*Q*TF;eZ=D*math.exp(D)
        for name,probe in block['probes'].items():
            meta=ref['probe_metadata'][name];VG=meta['variance'];TG=meta['tail_variance_bound']
            eL=Q*math.sqrt(TF*TG);L=abs(meta['mean'])+Q*math.sqrt(VF*VG)
            eA=eZ*L+(1+eZ)*eL
            eB=eZ*(L*L+VG)+(1+eZ)*(2*L*eL+eL*eL+TG)
            bounds={'matrix_first':eA,'matrix_second':eB}
            for state,m in probe['states'].items():
                denominator=m['probability']-eZ
                if denominator<=0:raise ValueError('Truncation cannot bound branch normalization')
                bounds[state+'/mean']=(eA+abs(m['mean'])*eZ)/denominator
                bounds[state+'/second']=(eB+abs(m['second_moment'])*eZ)/denominator
            result[section+'/'+name]=bounds;maximum=max(maximum,*bounds.values())
    covariance=max(v for row in ref['probe_covariance_tail_bounds'].values() for v in row.values())
    return {'maximum':max(maximum,covariance),'bounds':result,'covariance_maximum':covariance,
      'scope':'Conservative Gaussian bounds evaluated in floating point, not interval certified.'}

def normalization_checks(d):
    maximum=0.;minimum_single=1.;minimum_double=1.
    for section,block in sections(d).items():
        gram=decode_matrix(block['gram'])
        maximum=max(maximum,float(np.max(abs(np.diag(gram)-1))),float(np.max(abs(gram-gram.conj().T))))
        for probe in block['probes'].values():
            states=probe['states'];keys=('plus','minus') if section.startswith('single') else tuple(states)
            maximum=max(maximum,abs(sum(states[s]['probability'] for s in keys)-1))
            if section.startswith('single'):minimum_single=min(minimum_single,*(states[s]['probability'] for s in keys))
            else:minimum_double=min(minimum_double,*(states[s]['probability'] for s in keys))
            for m in states.values():
                if not all(math.isfinite(v) for v in m.values()):raise ValueError('Nonfinite state moment')
                if m['variance'] < -1e-9:raise ValueError('Negative variance')
    return {'maximum_error':float(maximum),'minimum_single_probability':float(minimum_single),
      'minimum_double_probability':float(minimum_double)}

def numeric_flatten(value,prefix=''):
    out={}
    if isinstance(value,dict):
        for k,v in value.items():out.update(numeric_flatten(v,prefix+'/'+str(k)))
    elif isinstance(value,(list,tuple)):
        for i,v in enumerate(value):out.update(numeric_flatten(v,prefix+'/'+str(i)))
    elif isinstance(value,(int,float,complex,np.number)) and not isinstance(value,bool):
        if not np.isfinite(value):raise ValueError('Nonfinite diagnostic value at '+prefix)
        out[prefix]=value
    return out

def numeric_difference(left,right):
    a,b=numeric_flatten(left),numeric_flatten(right)
    if set(a)!=set(b):raise ValueError('Refinement payload schemas differ')
    errors={k:float(abs(a[k]-b[k])) for k in a}
    worst=max(errors,key=errors.get) if errors else None
    return {'maximum':errors[worst] if worst else 0.,'worst':worst,'errors':errors}

def response_payload(d,diagnostic_key):
    diag=d[diagnostic_key]
    return {**{k:d[k] for k in ('means','covariance','first_moment','responses')},
      'centroid_spectrum':diag['centroid_spectrum'],
      'whitened_responses':diag['whitened_responses'],
      'whitened_cos_eigenvalues':diag['whitened_cos_eigenvalues'],
      'galerkin_residuals':diag['galerkin_residuals']}

def lattice_payload(d):
    return {**response_payload(d,'diagnostics'),
      **{k:d[k] for k in ('single_current_rules','coincident_double_current_rules',
        'twice_kick_current_rule','delayed_double_current_rules','instant_kick_scan')},
      'instant_observables':moment_values(d['instant_observables'])}

def relative_whitened_error(value,target):
    A=np.asarray(value,dtype=float);B=np.asarray(target,dtype=float)
    w,U=np.linalg.eigh((B+B.T)/2)
    if w[0]<=0:raise ValueError('Reference static metric is not positive definite')
    inverse=(U*(1/np.sqrt(w)))@U.T;D=inverse@(A-B)@inverse
    return {'operator_norm':float(np.linalg.norm(D,2)),'whitened_difference':D.tolist()}

def pair(z):return [float(np.real(z)),float(np.imag(z))]

def per_point(d,cfg):
    ref=d['refined_reference'];D=d['diagnostics'];R=ref['primary_diagnostics']
    C=np.asarray(d['covariance']);Cc=np.asarray(ref['covariance']);K=np.asarray(d['first_moment']);Kc=np.asarray(ref['first_moment'])
    covariance=relative_whitened_error(C[:3,:3],Cc[:3,:3]);first=relative_whitened_error(K[:3,:3],Kc[:3,:3])
    la=np.asarray(D['centroid_spectrum']);lc=np.asarray(R['centroid_spectrum'])
    if np.any(lc<=0):raise ValueError('Nonpositive reference centroid')
    centroid=np.abs(la-lc)/lc;cosine={};autos={};shape_max=0.;closure=0.
    Va=np.diag(C);Vc=np.diag(Cc)
    if np.any(Va<=0) or np.any(Vc<=0):raise ValueError('Nonpositive source variance')
    for t,encoded in d['responses'].items():
        A=decode_matrix(encoded);B=decode_matrix(ref['responses'][t]);ra=np.diag(A)/Va;rc=np.diag(B)/Vc
        cosine[t]=np.abs(np.asarray(D['whitened_cos_eigenvalues'][t])-R['whitened_cos_eigenvalues'][t]).tolist()
        autos[t]={}
        for i,name in enumerate(d['probe_names']):
            difference=A[i,i]-B[i,i];static1=(Va[i]-Vc[i])*rc[i];shape1=Va[i]*(ra[i]-rc[i])
            static2=(Va[i]-Vc[i])*ra[i];shape2=Vc[i]*(ra[i]-rc[i])
            error=float(abs(ra[i]-rc[i]));shape_max=max(shape_max,error)
            closure=max(closure,float(abs(difference-static1-shape1)),float(abs(difference-static2-shape2)))
            autos[t][name]={'lattice_variance':float(Va[i]),'continuum_variance':float(Vc[i]),
              'variance_ratio':float(Va[i]/Vc[i]),'mean_offset':float(d['means'][i]-ref['means'][i]),
              'lattice_normalized_response':pair(ra[i]),'continuum_normalized_response':pair(rc[i]),
              'normalized_response_difference':pair(ra[i]-rc[i]),'normalized_response_absolute_error':error,
              'raw_response_difference':pair(difference),
              'continuum_shape_order':{'static_weight':pair(static1),'normalized_shape':pair(shape1)},
              'lattice_shape_order':{'static_weight':pair(static2),'normalized_shape':pair(shape2)}}
    Sc=np.asarray(ref['fsum'])[:3,:3];S0=np.asarray(d['single_current_rules']['vacuum']['fsum'])[:3,:3]
    def current_records(group):
        records={}
        for state,row in group.items():
            S=np.asarray(row['fsum'])[:3,:3]
            absolute=relative_whitened_error(S,Sc);disturbance=relative_whitened_error(S,S0)
            records[state]={'relative_to_reference':absolute['operator_norm'],
              'disturbance_from_vacuum':disturbance['operator_norm'],'primary_fsum':S.tolist(),
              'reference_whitened_difference':absolute['whitened_difference'],
              'vacuum_whitened_difference':disturbance['whitened_difference']}
            for key in ('probability','maximum_link_phase'):
                if key in row:records[state][key]=row[key]
        return records
    currents={k:current_records(d[k]) for k in ('single_current_rules','coincident_double_current_rules',
      'twice_kick_current_rule','delayed_double_current_rules')}
    scan=current_records(d['instant_kick_scan']);instant=moment_errors(d['instant_observables'],d['instant_refined_reference'])
    return {'static_covariance_relative':covariance['operator_norm'],'static_first_moment_relative':first['operator_norm'],
      'static_covariance':covariance,'static_first_moment':first,
      'max_normalized_auto_response_error':shape_max,'auto_response_decompositions':autos,'decomposition_closure_error':closure,
      'centroid_relative_error':float(np.max(centroid)),'centroid_relative_errors':centroid.tolist(),
      'lattice_centroids':la.tolist(),'continuum_centroids':lc.tolist(),
      'cosine_eigen_absolute_error':max(v for row in cosine.values() for v in row),'cosine_eigen_errors':cosine,
      'current_rules':currents,'instant_kick_scan':scan,'reference_primary_fsum':Sc.tolist(),
      'instantaneous_moments':{'errors':instant,'maximum':max(instant.values()),'worst':max(instant,key=instant.get)},
      'lattice_galerkin_residuals':{t:r['whitened_operator_norm'] for t,r in D['galerkin_residuals'].items()},
      'continuum_galerkin_residuals':{t:r['whitened_operator_norm'] for t,r in R['galerkin_residuals'].items()},
      'max_lattice_galerkin_residual':max(r['whitened_operator_norm'] for r in D['galerkin_residuals'].values()),
      'max_continuum_galerkin_residual':max(r['whitened_operator_norm'] for r in R['galerkin_residuals'].values())}

def analyze(out):
    out=Path(out);cfg=card();g=cfg['gates']
    absent=set(GATE_DEFAULTS)-set(g)
    if absent:raise ValueError('Missing frozen gate keys: '+','.join(sorted(absent)))
    labels=[f'N{n}' for n in cfg['N']]+['N22_refined'];missing=[x for x in labels if not (out/(x+'.json')).exists()]
    if not (out/'prechecks.json').exists():missing.append('prechecks')
    if missing:return {'complete':False,'pass':False,'numerically_valid':False,'verdict':'incomplete',
      'missing':missing,'massive_status':'withheld','theory_status':'not established'}
    points={x:json.loads((out/(x+'.json')).read_text()) for x in labels};pre=json.loads((out/'prechecks.json').read_text())
    numeric_flatten(points)
    inputs=json.loads((HERE/'INPUTS.json').read_text())
    metrics={x:per_point(d,cfg) for x,d in points.items()}
    normalization={x:{'lattice':normalization_checks(d['instant_observables']),
      'continuum':normalization_checks(d['instant_refined_reference'])} for x,d in points.items()}
    tails={};doubling={};instant_doubling={}
    for x,d in points.items():
        ref=d['refined_reference'];instant=instant_tail_control(d['instant_refined_reference'],cfg)
        E=ref['primary_covariance_tail_norm'];eigen=ref['primary_diagnostics']['gram_eigenvalues'];cmin=eigen[0]
        generator_tail=E*float(np.linalg.norm(np.asarray(ref['first_moment'])[:3,:3],2))/(cmin*cmin)
        minimum_centroid=min(ref['primary_diagnostics']['centroid_spectrum'])
        if minimum_centroid<=generator_tail:raise ValueError('Reference centroid tail overlaps zero')
        covariance_relative_tail=E*float(np.linalg.norm(np.asarray(d['covariance'])[:3,:3],2))/(cmin*cmin)
        centroid_relative_tail=(1+metrics[x]['centroid_relative_error'])*generator_tail/(minimum_centroid-generator_tail)
        galerkin_tail=ref['primary_whitened_response_tail_bound']+max(abs(t) for t in cfg['response_times'])*generator_tail
        tails[x]={'instantaneous':instant,'whitened_response':ref['primary_whitened_response_tail_bound'],
          'static_covariance_relative':covariance_relative_tail,'projected_generator':generator_tail,
          'centroid_relative':centroid_relative_tail,'galerkin_residual':galerkin_tail,
          'normalized_auto_response':max(ref['normalized_auto_response_tail_bounds']),
          'raw_response_entry':float(np.max(ref['entry_tail_bounds']))}
        tails[x]['maximum']=max(instant['maximum'],tails[x]['whitened_response'],tails[x]['normalized_auto_response'],
          tails[x]['raw_response_entry'],covariance_relative_tail,centroid_relative_tail,galerkin_tail)
        doubling[x]=numeric_difference(response_payload(d['reference'],'primary_diagnostics'),response_payload(ref,'primary_diagnostics'))
        ie=moment_errors(d['instant_reference'],d['instant_refined_reference'])
        instant_doubling[x]={'maximum':max(ie.values()),'worst':max(ie,key=ie.get),'errors':ie}
    refinement=numeric_difference(lattice_payload(points['N22']),lattice_payload(points['N22_refined']))
    numerical_fields=('input_norm_error','ground_residual','max_static_imaginary','max_free_norm_error',
      'sumrule_identity_error','kick_cosine_identity_error','twice_kick_cosine_identity_error')
    rule_probabilities={}
    for x,d in points.items():
        single=d['single_current_rules'];coincident=d['coincident_double_current_rules'];delayed=d['delayed_double_current_rules']
        rule_probabilities[x]={'maximum_error':max(abs(single['vacuum']['probability']-1),abs(single['kick']['probability']-1),
          abs(single['plus']['probability']+single['minus']['probability']-1),
          abs(sum(r['probability'] for r in coincident.values())-1),abs(sum(r['probability'] for r in delayed.values())-1),
          abs(d['twice_kick_current_rule']['double_kick']['probability']-1)),
          'minimum_single':min(single[s]['probability'] for s in ('plus','minus')),
          'minimum_double':min(r['probability'] for block in (coincident,delayed) for r in block.values())}
    quality={'source_and_reference_prechecks':bool(pre['passed']),
      'registered_points':all(d['N']==int(x[1:].split('_')[0]) and d['dimension']==math.comb(d['N'],d['N']//2)
        and abs(d['a']-cfg['length']/(d['N']+1))<1e-14 and d['refined']==x.endswith('_refined')
        and d['input_sha256']==inputs[x]['sha256'] for x,d in points.items()),
      'independent_source_weights':all(float(np.max(abs(np.asarray(d['probe_weights'])-ref['probe_weights'])))<=g['numerical']
        and max(ref['probe_partition_errors'])<=g['numerical'] for d in points.values() for ref in (d['reference'],d['refined_reference'])),
      'numerical':all(max(d[k] for k in numerical_fields)<=g['numerical'] for d in points.values()),
      'normalization':all(r['maximum_error']<=g['numerical'] for block in normalization.values() for r in block.values())
        and all(r['maximum_error']<=g['numerical'] for r in rule_probabilities.values()),
      'conditioning':all(r['minimum_single_probability']>=g['single_probability_floor'] and r['minimum_double_probability']>=g['double_probability_floor'] for block in normalization.values() for r in block.values())
        and all(r['minimum_single']>=g['single_probability_floor'] and r['minimum_double']>=g['double_probability_floor'] for r in rule_probabilities.values()),
      'reference_imaginary':all(ref['max_imaginary_residual']<=g['numerical'] for d in points.values() for ref in (d['instant_reference'],d['instant_refined_reference'])),
      'reference_mode_doubling':max(r['maximum'] for r in doubling.values())<=g['reference_doubling'],
      'instant_reference_mode_doubling':max(r['maximum'] for r in instant_doubling.values())<=g['reference_doubling'],
      'reference_tail':max(r['maximum'] for r in tails.values())<=g['reference_tail'],
      'lattice_refinement':refinement['maximum']<=g['refinement'],
      'signed_decomposition_identity':max(m['decomposition_closure_error'] for m in metrics.values())<=g['numerical'],
      'response_contraction':all(d[k]['maximum_response_contraction_excess']<=g['numerical'] for p in points.values() for d,k in ((p,'diagnostics'),(p['reference'],'primary_diagnostics'),(p['refined_reference'],'primary_diagnostics'))),
      'moment_hermiticity':all(max(d[k]['gram_hermiticity_error'],d[k]['first_moment_hermiticity_error'])<=g['numerical'] for p in points.values() for d,k in ((p,'diagnostics'),(p['reference'],'primary_diagnostics'),(p['refined_reference'],'primary_diagnostics'))),
      'primary_gram_conditioning':all(d[k]['gram_eigenvalues'][0]>cfg['gram_minimum_eigenvalue'] and d[k]['gram_condition']<cfg['gram_condition_limit'] for p in points.values() for d,k in ((p,'diagnostics'),(p['reference'],'primary_diagnostics'),(p['refined_reference'],'primary_diagnostics')))}
    valid=all(quality.values());finest=('N22','N22_refined')
    def hypothesis(values,threshold,guard=None):
        guard=g['numerical_guard'] if guard is None else guard
        accepted=all(v+guard<=threshold for v in values.values())
        return {'threshold':threshold,'guard':guard,'finest_metrics':values,'guarded_within_allowance':accepted,
          'passed':bool(valid and accepted),'status':('supported_at_registered_points' if accepted else 'not_supported_at_registered_points') if valid else 'numerically_unresolved'}
    specs={'static_covariance':('static_covariance_relative','static_relative'),
      'static_first_moment':('static_first_moment_relative','static_relative'),
      'normalized_auto_response':('max_normalized_auto_response_error','normalized_response_absolute'),
      'centroid_spectrum':('centroid_relative_error','centroid_relative'),
      'cosine_spectrum':('cosine_eigen_absolute_error','cosine_eigen_absolute'),
      'lattice_source_span_predictor':('max_lattice_galerkin_residual','galerkin_residual'),
      'continuum_source_span_predictor':('max_continuum_galerkin_residual','galerkin_residual')}
    hypotheses={name:hypothesis({x:metrics[x][field] for x in finest},g[gate]) for name,(field,gate) in specs.items()}
    hypotheses['instantaneous_moments']=hypothesis({x:metrics[x]['instantaneous_moments']['maximum'] for x in finest},g['instant_moment_absolute'],g['instant_moment_reserve'])
    for group in ('single_current_rules','coincident_double_current_rules','twice_kick_current_rule','delayed_double_current_rules'):
        hypotheses[group]=hypothesis({x:max(r['relative_to_reference'] for r in metrics[x]['current_rules'][group].values()) for x in finest},g['current_relative'])
    scan={}
    for coupling in points['N22']['instant_kick_scan']:
        disturbance=hypothesis({x:metrics[x]['instant_kick_scan'][coupling]['disturbance_from_vacuum'] for x in finest},g['current_relative'])
        absolute=hypothesis({x:metrics[x]['instant_kick_scan'][coupling]['relative_to_reference'] for x in finest},g['current_relative'])
        scan[coupling]={'disturbance':disturbance,'absolute_current':absolute,'both_supported':disturbance['passed'] and absolute['passed']}
    trend={x:{k:m[k] for k in ('static_covariance_relative','static_first_moment_relative',
      'max_normalized_auto_response_error','centroid_relative_error','cosine_eigen_absolute_error',
      'max_lattice_galerkin_residual','max_continuum_galerkin_residual')} for x,m in metrics.items()}
    return {'complete':True,'pass':valid,'numerically_valid':valid,
      'verdict':'numerically_valid_diagnostic_hypotheses_reported_separately' if valid else 'numerical_validity_unresolved',
      'quality_gates':quality,'hypotheses':hypotheses,'instant_kick_scan_hypotheses':scan,
      'metrics':metrics,'reference_tails':tails,'normalization':normalization,'current_probability_checks':rule_probabilities,
      'reference_mode_doubling':doubling,'instant_reference_mode_doubling':instant_doubling,'refinement':refinement,'trend':trend,
      'massive_status':'withheld','N22_historical_verdict_unchanged':True,'theory_status':'not established',
      'scope':cfg['scope'],'interpretation':'Pass denotes numerical validity only. Static normalization, source-span invariants, current rules and projected evolution are separate finite diagnostics. No unique attribution to Hamiltonian evolution follows from a failed normalized response, and no theory or empirical-validation pass is defined.'}
