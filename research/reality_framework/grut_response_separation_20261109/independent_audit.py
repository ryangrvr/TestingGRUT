"""Independent arithmetic audit of archived R1 outputs; no evolution imports.

Run only after the prospective protocol has frozen and production has finished.
This does not re-run evolution, audit binary checkpoints, or certify the continuum.
"""
import argparse
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

AREA = Path('research/reality_framework/grut_response_separation_20261109')


def array(x):
    return np.asarray(x, dtype=float)


def complex_array(x):
    a = array(x)
    return a[..., 0] + 1j * a[..., 1]


def hermitian(x):
    return (x + x.conj().T) / 2


def roots(x):
    values, vectors = np.linalg.eigh(hermitian(x))
    if values[0] <= 0:
        raise ValueError('Positive covariance/sum-rule metric required')
    return ((vectors * np.sqrt(values)) @ vectors.conj().T,
            (vectors / np.sqrt(values)) @ vectors.conj().T)


def opnorm(x):
    return float(np.linalg.svd(x, compute_uv=False)[0])


def whitened_difference(metric, difference):
    _, w = roots(metric)
    return opnorm(w @ difference @ w)


def encoded(matrix):
    return np.stack((matrix.real,matrix.imag),axis=-1).tolist()


def independent_diagnostics(covariance,first_moment,responses):
    c = array(covariance)[:3,:3]
    k = array(first_moment)[:3,:3]
    q,w = roots(c)
    a = hermitian(w@k@w)
    eigen = np.linalg.eigvalsh(c)
    result = {'gram_eigenvalues':eigen.tolist(),'gram_condition':float(eigen[-1]/eigen[0]),
        'gram_hermiticity_error':float(np.max(abs(c-c.T))),
        'first_moment_hermiticity_error':float(np.max(abs(k-k.T))),
        'projected_generator':encoded(a),'centroid_spectrum':np.linalg.eigvalsh(a).tolist(),
        'whitened_responses':{},'whitened_cos_eigenvalues':{},'galerkin_responses':{},
        'whitened_galerkin_responses':{},'galerkin_residuals':{}}
    excess = 0.
    for key,value in responses.items():
        ct = complex_array(value)[:3,:3]
        normalized = w@ct@w
        white_prediction = expm(-1j*float(key)*a)
        prediction = q@white_prediction@q
        delta = normalized-white_prediction
        result['whitened_responses'][key] = encoded(normalized)
        result['whitened_cos_eigenvalues'][key] = np.linalg.eigvalsh(hermitian(normalized)).tolist()
        result['galerkin_responses'][key] = encoded(prediction)
        result['whitened_galerkin_responses'][key] = encoded(white_prediction)
        result['galerkin_residuals'][key] = {'whitened_operator_norm':opnorm(delta),
            'whitened_matrix':encoded(delta),'normalized_diagonal_absolute':(abs(np.diag(ct-prediction))/np.diag(c)).tolist()}
        excess = max(excess,opnorm(normalized)-1)
    result['maximum_response_contraction_excess'] = excess
    return result


def observation_blocks(record):
    return {**{'single/'+key:value for key,value in record['single'].items()},'double':record['double']}


def observation_values(record):
    values = {}
    for label,block in observation_blocks(record).items():
        for name,probe in block['probes'].items():
            for state,row in probe['states'].items():
                for field in ('mean','second_moment'):
                    values[f'{label}/{name}/{state}/{field}'] = row[field]
            for field in ('first_matrix','second_matrix'):
                for i,row in enumerate(complex_array(probe[field])):
                    for j,value in enumerate(row):
                        values[f'{label}/{name}/{field}/{i}/{j}'] = value
    for name,row in record['probe_covariance'].items():
        for other,value in row.items():
            values[f'vacuum_covariance/{name}/{other}'] = value
    return values


def observation_errors(left,right):
    a,b = observation_values(left),observation_values(right)
    if a.keys() != b.keys():
        raise ValueError('Observable dictionaries differ')
    return {key:float(abs(value-b[key])) for key,value in a.items()}


def recompute_metrics(point,cfg):
    ref = point['refined_reference']
    c,k = array(point['covariance']),array(point['first_moment'])
    cc,kc = array(ref['covariance']),array(ref['first_moment'])
    d = independent_diagnostics(c,k,point['responses'])
    r = independent_diagnostics(cc,kc,ref['responses'])
    def static(value,target):
        _,w = roots(target)
        delta = w@(value-target)@w
        return {'operator_norm':opnorm(delta),'whitened_difference':delta.tolist()}
    static_c,static_k = static(c[:3,:3],cc[:3,:3]),static(k[:3,:3],kc[:3,:3])
    va,vc = np.diag(c),np.diag(cc)
    autos = {};cosines = {};closure = 0.;maxshape = 0.
    def pair(value):return [float(value.real),float(value.imag)]
    for key,raw in point['responses'].items():
        left,right = complex_array(raw),complex_array(ref['responses'][key])
        ra,rc = np.diag(left)/va,np.diag(right)/vc
        autos[key] = {}
        cosines[key] = abs(array(d['whitened_cos_eigenvalues'][key])-array(r['whitened_cos_eigenvalues'][key])).tolist()
        for i,name in enumerate(point['probe_names']):
            delta = left[i,i]-right[i,i]
            static1,shape1 = (va[i]-vc[i])*rc[i],va[i]*(ra[i]-rc[i])
            static2,shape2 = (va[i]-vc[i])*ra[i],vc[i]*(ra[i]-rc[i])
            shape_error = float(abs(ra[i]-rc[i]))
            maxshape = max(maxshape,shape_error)
            closure = max(closure,float(abs(delta-static1-shape1)),float(abs(delta-static2-shape2)))
            autos[key][name] = {'lattice_variance':float(va[i]),'continuum_variance':float(vc[i]),
                'variance_ratio':float(va[i]/vc[i]),'mean_offset':float(point['means'][i]-ref['means'][i]),
                'lattice_normalized_response':pair(ra[i]),'continuum_normalized_response':pair(rc[i]),
                'normalized_response_difference':pair(ra[i]-rc[i]),'normalized_response_absolute_error':shape_error,
                'raw_response_difference':pair(delta),
                'continuum_shape_order':{'static_weight':pair(static1),'normalized_shape':pair(shape1)},
                'lattice_shape_order':{'static_weight':pair(static2),'normalized_shape':pair(shape2)}}
    sc,s0 = array(ref['fsum'])[:3,:3],array(point['single_current_rules']['vacuum']['fsum'])[:3,:3]
    def current_group(group):
        records = {}
        for label,row in group.items():
            s = array(row['fsum'])[:3,:3]
            records[label] = {'relative_to_reference':whitened_difference(sc,s-sc),
                'disturbance_from_vacuum':whitened_difference(s0,s-s0),'primary_fsum':s.tolist()}
            for key in ('probability','maximum_link_phase'):
                if key in row:records[label][key] = row[key]
        return records
    current_keys = ('single_current_rules','coincident_double_current_rules','twice_kick_current_rule','delayed_double_current_rules')
    currents = {key:current_group(point[key]) for key in current_keys}
    instant = observation_errors(point['instant_observables'],point['instant_refined_reference'])
    la,lc = array(d['centroid_spectrum']),array(r['centroid_spectrum'])
    centroids = abs(la-lc)/lc
    lattice_gal = {key:row['whitened_operator_norm'] for key,row in d['galerkin_residuals'].items()}
    continuum_gal = {key:row['whitened_operator_norm'] for key,row in r['galerkin_residuals'].items()}
    return {'static_covariance_relative':static_c['operator_norm'],'static_first_moment_relative':static_k['operator_norm'],
        'static_covariance':static_c,'static_first_moment':static_k,
        'max_normalized_auto_response_error':maxshape,'auto_response_decompositions':autos,'decomposition_closure_error':closure,
        'centroid_relative_error':float(max(centroids)),'centroid_relative_errors':centroids.tolist(),
        'lattice_centroids':la.tolist(),'continuum_centroids':lc.tolist(),
        'cosine_eigen_absolute_error':max(value for row in cosines.values() for value in row),'cosine_eigen_errors':cosines,
        'current_rules':currents,'instant_kick_scan':current_group(point['instant_kick_scan']),'reference_primary_fsum':sc.tolist(),
        'instantaneous_moments':{'errors':instant,'maximum':max(instant.values()),'worst':max(instant,key=instant.get)},
        'lattice_galerkin_residuals':lattice_gal,'continuum_galerkin_residuals':continuum_gal,
        'max_lattice_galerkin_residual':max(lattice_gal.values()),'max_continuum_galerkin_residual':max(continuum_gal.values())}


def flatten(value,prefix=''):
    records = {}
    if isinstance(value,dict):
        for key,item in value.items():records.update(flatten(item,prefix+'/'+str(key)))
    elif isinstance(value,(tuple,list)):
        for key,item in enumerate(value):records.update(flatten(item,prefix+'/'+str(key)))
    elif isinstance(value,(int,float,complex,np.number)) and not isinstance(value,bool):
        if not np.isfinite(value):raise ValueError('Nonfinite number: '+prefix)
        records[prefix] = value
    return records


def payload_difference(left,right):
    a,b = flatten(left),flatten(right)
    if a.keys() != b.keys():raise ValueError('Refinement schemas differ')
    errors = {key:float(abs(value-b[key])) for key,value in a.items()}
    worst = max(errors,key=errors.get)
    return {'maximum':errors[worst],'worst':worst,'errors':errors}


def response_values(record,diagnostic_key):
    diagnostics = record[diagnostic_key]
    return {**{key:record[key] for key in ('means','covariance','first_moment','responses')},
        **{key:diagnostics[key] for key in ('centroid_spectrum','whitened_responses',
            'whitened_cos_eigenvalues','galerkin_residuals')}}


def lattice_values(record):
    return {**response_values(record,'diagnostics'),
        **{key:record[key] for key in ('single_current_rules','coincident_double_current_rules',
            'twice_kick_current_rule','delayed_double_current_rules','instant_kick_scan')},
        'instant_observables':observation_values(record['instant_observables'])}


def normalization(record):
    errors = [0.];single = [1.];double = [1.]
    for label,block in observation_blocks(record).items():
        gram = complex_array(block['gram'])
        errors.extend([float(np.max(abs(np.diag(gram)-1))),float(np.max(abs(gram-gram.conj().T)))])
        for row in block['probes'].values():
            labels = ('plus','minus') if label.startswith('single/') else tuple(row['states'])
            probabilities = [row['states'][key]['probability'] for key in labels]
            errors.append(abs(sum(probabilities)-1))
            (single if label.startswith('single/') else double).extend(probabilities)
    return {'maximum_error':max(errors),'minimum_single_probability':min(single),
        'minimum_double_probability':min(double)}


def instant_tails(record,cfg):
    detector_variance = record['detector']['variance']
    detector_tail = record['detector']['tail_variance_bound']
    bounds = {}
    for label,block in observation_blocks(record).items():
        charge = cfg['coupling']*(4 if label=='double' else 2)
        exponent = charge**2*detector_tail/2
        gram_tail = exponent*math.exp(exponent)
        for name,probe in block['probes'].items():
            metadata = record['probe_metadata'][name]
            variance,tail = metadata['variance'],metadata['tail_variance_bound']
            shift_tail = charge*math.sqrt(detector_tail*tail)
            max_shift = abs(metadata['mean'])+charge*math.sqrt(detector_variance*variance)
            first_tail = gram_tail*max_shift+(1+gram_tail)*shift_tail
            second_tail = gram_tail*(max_shift**2+variance)+(1+gram_tail)*(2*max_shift*shift_tail+shift_tail**2+tail)
            row = {'matrix_first':first_tail,'matrix_second':second_tail}
            for state,value in probe['states'].items():
                denominator = value['probability']-gram_tail
                if denominator<=0:raise ValueError('Uncontrolled conditional tail')
                row[state+'/mean'] = (first_tail+abs(value['mean'])*gram_tail)/denominator
                row[state+'/second'] = (second_tail+abs(value['second_moment'])*gram_tail)/denominator
            bounds[label+'/'+name] = row
    covariance = max(value for row in record['probe_covariance_tail_bounds'].values() for value in row.values())
    return {'maximum':max(covariance,*(value for row in bounds.values() for value in row.values())),
        'bounds':bounds,'covariance_maximum':covariance}


def recompute_assessment(points,cfg,manifest,pre):
    g = cfg['gates'];metrics = {label:recompute_metrics(point,cfg) for label,point in points.items()}
    norm = {};probabilities = {};tails = {};doubling = {};instant_doubling = {}
    for label,point in points.items():
        norm[label] = {'lattice':normalization(point['instant_observables']),
            'continuum':normalization(point['instant_refined_reference'])}
        single,coincident,delayed = (point[key] for key in ('single_current_rules',
            'coincident_double_current_rules','delayed_double_current_rules'))
        probabilities[label] = {'maximum_error':max(abs(single['vacuum']['probability']-1),
            abs(single['kick']['probability']-1),abs(single['plus']['probability']+single['minus']['probability']-1),
            abs(sum(row['probability'] for row in coincident.values())-1),
            abs(sum(row['probability'] for row in delayed.values())-1),
            abs(point['twice_kick_current_rule']['double_kick']['probability']-1)),
            'minimum_single':min(single[key]['probability'] for key in ('plus','minus')),
            'minimum_double':min(row['probability'] for group in (coincident,delayed) for row in group.values())}
        ref = point['refined_reference']
        # Recompute source variation and covariance tail from archived cell weights.
        heights = array(ref['probe_weights'])/ref['a']
        variation = abs(heights[:,0])+np.sum(abs(np.diff(heights,axis=1)),axis=1)+abs(heights[:,-1])
        variance_tails = variation**2*cfg['length']**2/(2*np.pi**4*ref['modes']**2)
        entry_tail = np.sqrt(np.outer(variance_tails,variance_tails))
        E = opnorm(entry_tail[:3,:3])
        eigen = np.linalg.eigvalsh(array(ref['covariance'])[:3,:3]);cmin,cmax = eigen[0],eigen[-1]
        white_tail = E/cmin+E*cmax/(cmin*cmin)
        generator_tail = E*opnorm(array(ref['first_moment'])[:3,:3])/(cmin*cmin)
        covariance_tail = E*opnorm(array(point['covariance'])[:3,:3])/(cmin*cmin)
        centroid_min = min(independent_diagnostics(ref['covariance'],ref['first_moment'],ref['responses'])['centroid_spectrum'])
        centroid_tail = (1+metrics[label]['centroid_relative_error'])*generator_tail/(centroid_min-generator_tail)
        gal_tail = white_tail+max(abs(t) for t in cfg['response_times'])*generator_tail
        instantaneous = instant_tails(point['instant_refined_reference'],cfg)
        tails[label] = {'instantaneous':instantaneous,'whitened_response':white_tail,
            'static_covariance_relative':covariance_tail,'projected_generator':generator_tail,
            'centroid_relative':centroid_tail,'galerkin_residual':gal_tail,
            'normalized_auto_response':float(max(2*variance_tails/np.diag(array(ref['covariance'])))),
            'raw_response_entry':float(np.max(entry_tail))}
        tails[label]['maximum'] = max(instantaneous['maximum'],*(value for key,value in tails[label].items()
            if key not in ('instantaneous','projected_generator')))
        doubling[label] = payload_difference(response_values(point['reference'],'primary_diagnostics'),
            response_values(ref,'primary_diagnostics'))
        errors = observation_errors(point['instant_reference'],point['instant_refined_reference'])
        worst = max(errors,key=errors.get)
        instant_doubling[label] = {'maximum':errors[worst],'worst':worst,'errors':errors}
    refinement = payload_difference(lattice_values(points['N22']),lattice_values(points['N22_refined']))
    numerical_fields = ('input_norm_error','ground_residual','max_static_imaginary','max_free_norm_error',
        'sumrule_identity_error','kick_cosine_identity_error','twice_kick_cosine_identity_error')
    diagnostics = [independent_diagnostics(side['covariance'],side['first_moment'],side['responses'])
        for point in points.values() for side in (point,point['reference'],point['refined_reference'])]
    quality = {
        'source_and_reference_prechecks':bool(pre['passed']),
        'registered_points':all(point['N']==int(label[1:].split('_')[0]) and point['dimension']==math.comb(point['N'],point['N']//2)
            and abs(point['a']-cfg['length']/(point['N']+1))<1e-14 and point['refined']==label.endswith('_refined')
            and point['input_sha256']==manifest[label]['sha256'] for label,point in points.items()),
        'independent_source_weights':all(np.max(abs(array(point['probe_weights'])-array(ref['probe_weights'])))<=g['numerical']
            and max(ref['probe_partition_errors'])<=g['numerical'] for point in points.values() for ref in (point['reference'],point['refined_reference'])),
        'numerical':all(max(point[key] for key in numerical_fields)<=g['numerical'] for point in points.values()),
        'normalization':all(row['maximum_error']<=g['numerical'] for group in norm.values() for row in group.values())
            and all(row['maximum_error']<=g['numerical'] for row in probabilities.values()),
        'conditioning':all(row['minimum_single_probability']>=g['single_probability_floor'] and row['minimum_double_probability']>=g['double_probability_floor']
            for group in norm.values() for row in group.values()) and all(row['minimum_single']>=g['single_probability_floor']
            and row['minimum_double']>=g['double_probability_floor'] for row in probabilities.values()),
        'reference_imaginary':all(ref['max_imaginary_residual']<=g['numerical'] for point in points.values()
            for ref in (point['instant_reference'],point['instant_refined_reference'])),
        'reference_mode_doubling':max(row['maximum'] for row in doubling.values())<=g['reference_doubling'],
        'instant_reference_mode_doubling':max(row['maximum'] for row in instant_doubling.values())<=g['reference_doubling'],
        'reference_tail':max(row['maximum'] for row in tails.values())<=g['reference_tail'],
        'lattice_refinement':refinement['maximum']<=g['refinement'],
        'signed_decomposition_identity':max(row['decomposition_closure_error'] for row in metrics.values())<=g['numerical'],
        'response_contraction':all(row['maximum_response_contraction_excess']<=g['numerical'] for row in diagnostics),
        'moment_hermiticity':all(max(row['gram_hermiticity_error'],row['first_moment_hermiticity_error'])<=g['numerical'] for row in diagnostics),
        'primary_gram_conditioning':all(row['gram_eigenvalues'][0]>cfg['gram_minimum_eigenvalue']
            and row['gram_condition']<cfg['gram_condition_limit'] for row in diagnostics)}
    valid = all(quality.values());finest = ('N22','N22_refined')
    def hypothesis(values,threshold,guard=None):
        guard = g['numerical_guard'] if guard is None else guard
        accepted = all(value+guard<=threshold for value in values.values())
        return {'threshold':threshold,'guard':guard,'finest_metrics':values,
            'guarded_within_allowance':accepted,'passed':bool(valid and accepted),
            'status':('supported_at_registered_points' if accepted else 'not_supported_at_registered_points') if valid else 'numerically_unresolved'}
    specs = {'static_covariance':('static_covariance_relative','static_relative'),
        'static_first_moment':('static_first_moment_relative','static_relative'),
        'normalized_auto_response':('max_normalized_auto_response_error','normalized_response_absolute'),
        'centroid_spectrum':('centroid_relative_error','centroid_relative'),
        'cosine_spectrum':('cosine_eigen_absolute_error','cosine_eigen_absolute'),
        'lattice_source_span_predictor':('max_lattice_galerkin_residual','galerkin_residual'),
        'continuum_source_span_predictor':('max_continuum_galerkin_residual','galerkin_residual')}
    hypotheses = {label:hypothesis({key:metrics[key][field] for key in finest},g[gate])
        for label,(field,gate) in specs.items()}
    hypotheses['instantaneous_moments'] = hypothesis({key:metrics[key]['instantaneous_moments']['maximum'] for key in finest},
        g['instant_moment_absolute'],g['instant_moment_reserve'])
    for group in ('single_current_rules','coincident_double_current_rules','twice_kick_current_rule','delayed_double_current_rules'):
        hypotheses[group] = hypothesis({key:max(row['relative_to_reference'] for row in metrics[key]['current_rules'][group].values())
            for key in finest},g['current_relative'])
    scan = {}
    for coupling in points['N22']['instant_kick_scan']:
        disturbance = hypothesis({key:metrics[key]['instant_kick_scan'][coupling]['disturbance_from_vacuum'] for key in finest},g['current_relative'])
        absolute = hypothesis({key:metrics[key]['instant_kick_scan'][coupling]['relative_to_reference'] for key in finest},g['current_relative'])
        scan[coupling] = {'disturbance':disturbance,'absolute_current':absolute,'both_supported':disturbance['passed'] and absolute['passed']}
    trend = {label:{key:row[key] for key in ('static_covariance_relative','static_first_moment_relative',
        'max_normalized_auto_response_error','centroid_relative_error','cosine_eigen_absolute_error',
        'max_lattice_galerkin_residual','max_continuum_galerkin_residual')} for label,row in metrics.items()}
    return {'complete':True,'pass':valid,'numerically_valid':valid,
        'verdict':'numerically_valid_diagnostic_hypotheses_reported_separately' if valid else 'numerical_validity_unresolved',
        'quality_gates':quality,'hypotheses':hypotheses,'instant_kick_scan_hypotheses':scan,'metrics':metrics,
        'reference_tails':tails,'normalization':norm,'current_probability_checks':probabilities,
        'reference_mode_doubling':doubling,'instant_reference_mode_doubling':instant_doubling,
        'refinement':refinement,'trend':trend,'massive_status':'withheld',
        'N22_historical_verdict_unchanged':True,'theory_status':'not established'}


def spectral_summary(covariance, first_moment, responses):
    c = array(covariance)[:3, :3]
    k = array(first_moment)[:3, :3]
    _, w = roots(c)
    generator = hermitian(w @ k @ w)
    result = {
        'covariance_eigenvalues': np.linalg.eigvalsh(c).tolist(),
        'centroids': np.linalg.eigvalsh(generator).tolist(),
        'condition': float(np.linalg.cond(c)),
        'times': {},
    }
    for key, raw in responses.items():
        t = float(key)
        ct = complex_array(raw)[:3, :3]
        normalized = w @ ct @ w
        result['times'][key] = {
            'cosine_eigenvalues': np.linalg.eigvalsh(hermitian(normalized)).tolist(),
            'galerkin_error': opnorm(normalized - expm(-1j * generator * t)),
            'compression_norm': opnorm(normalized),
        }
    return result


def audit_point(point, cfg):
    errors = {}
    def record(label, value):
        errors[label] = float(np.max(np.abs(value)))

    weights = array(point['probe_weights'])
    detector_weights = array(point['detector_weights'])
    c = array(point['covariance'])
    k = array(point['first_moment'])
    record('zero_response_covariance', complex_array(point['responses']['0']) - c)
    record('covariance_symmetry', c-c.T)
    record('first_moment_symmetry', k-k.T)
    groups = ('single_current_rules', 'coincident_double_current_rules',
              'twice_kick_current_rule', 'delayed_double_current_rules')
    for group in groups:
        for label, state in point[group].items():
            hops = array(state['hopping_expectations'])
            reconstructed = -(weights * hops) @ weights.T
            record(f'fsum/{group}/{label}', reconstructed-array(state['fsum']))
    vacuum = point['single_current_rules']['vacuum']
    hops0 = array(vacuum['hopping_expectations'])
    s0 = array(vacuum['fsum'])
    record('ground_fsum_identity', s0-2*k)
    kick = point['single_current_rules']['kick']
    twice = point['twice_kick_current_rule']['double_kick']
    record('kick_cosine', array(kick['hopping_expectations'])-hops0*np.cos(cfg['coupling']*detector_weights))
    record('twice_kick_cosine', array(twice['hopping_expectations'])-hops0*np.cos(2*cfg['coupling']*detector_weights))
    for key, scan in point['instant_kick_scan'].items():
        hops = hops0*np.cos(float(key)*detector_weights)
        record(f'kick_scan/{key}', -(weights*hops)@weights.T-array(scan['fsum']))
        record(f'kick_phase/{key}', np.max(abs(float(key)*detector_weights))-scan['maximum_link_phase'])

    single_coefficients = {
        'vacuum': np.array([1,0]), 'kick': np.array([0,1]),
        'plus': np.array([1,1])/2, 'minus': np.array([1,-1])/2,
    }
    double_coefficients = {
        ''.join('+' if s == 1 else '-' for s in signs):
        np.array([1, signs[0], signs[1], signs[0]*signs[1]])/4
        for signs in itertools.product((1,-1), repeat=2)
    }

    def audit_moments(label, observations):
        for family, records, coefficients in (
                ('single', observations['single']['0'], single_coefficients),
                ('double', observations['double'], double_coefficients)):
            gram = complex_array(records['gram'])
            record(f'{label}/{family}/gram_hermitian', gram-gram.conj().T)
            for name, probe in records['probes'].items():
                first = complex_array(probe['first_matrix'])
                second = complex_array(probe['second_matrix'])
                record(f'{label}/{family}/{name}/first_hermitian', first-first.conj().T)
                record(f'{label}/{family}/{name}/second_hermitian', second-second.conj().T)
                for state_name, x in coefficients.items():
                    state = probe['states'][state_name]
                    probability = x.conj()@gram@x
                    mean = (x.conj()@first@x)/probability
                    second_moment = (x.conj()@second@x)/probability
                    record(f'{label}/{family}/{name}/{state_name}/probability', probability-state['probability'])
                    record(f'{label}/{family}/{name}/{state_name}/mean', mean-state['mean'])
                    record(f'{label}/{family}/{name}/{state_name}/second_moment', second_moment-state['second_moment'])
                    record(f'{label}/{family}/{name}/{state_name}/variance', second_moment-mean**2-state['variance'])
                selective = ('plus','minus') if family == 'single' else tuple(double_coefficients)
                record(f'{label}/{family}/{name}/probability_sum', sum(probe['states'][s]['probability'] for s in selective)-1)
                if family == 'double':
                    for quantity in ('probability','mean','second_moment','variance'):
                        record(f'{label}/coincident_commutation/{name}/{quantity}',
                               probe['states']['+-'][quantity]-probe['states']['-+'][quantity])

    for key in ('instant_observables','instant_reference','instant_refined_reference'):
        audit_moments(key, point[key])
    for label, state in point['coincident_double_current_rules'].items():
        for name, probe in point['instant_observables']['double']['probes'].items():
            record(f'coincident_probability/{label}/{name}', state['probability']-probe['states'][label]['probability'])

    a = point['a']
    sc = (weights @ weights.T)/(np.pi*a)
    scan_metrics = {}
    for key, scan in point['instant_kick_scan'].items():
        s = array(scan['fsum'])[:3,:3]
        scan_metrics[key] = {
            'disturbance': whitened_difference(s0[:3,:3],s-s0[:3,:3]),
            'absolute_continuum': whitened_difference(sc[:3,:3],s-sc[:3,:3]),
        }
    state_metrics = {}
    for group in groups:
        state_metrics[group] = {}
        for label, state in point[group].items():
            s = array(state['fsum'])
            state_metrics[group][label] = {
                'diagonal_disturbance': (np.diag(s)/np.diag(s0)-1).tolist(),
                'full_primary_disturbance': whitened_difference(s0[:3,:3],s[:3,:3]-s0[:3,:3]),
            }
    return {
        'N': point['N'], 'refined': point['refined'],
        'largest_arithmetic_error': max(errors.values()), 'arithmetic_errors': errors,
        'independent_spectral_summary': spectral_summary(c,k,point['responses']),
        'independent_scan_metrics': scan_metrics, 'independent_state_metrics': state_metrics,
        'scope': 'Archived arithmetic audit only; no evolution or binary-checkpoint audit.',
    }


def compare_tree(saved, computed, path='', issues=None, maximum=None):
    if issues is None:
        issues = []
    if maximum is None:
        maximum = [0.]
    if isinstance(computed, dict):
        for key, value in computed.items():
            if key not in saved:
                issues.append(path+'/'+key+'/missing')
            else:
                compare_tree(saved[key],value,path+'/'+key,issues,maximum)
    elif isinstance(computed, (list,tuple)):
        if len(saved) != len(computed):
            issues.append(path+'/length')
        else:
            for index, (s,c) in enumerate(zip(saved,computed)):
                compare_tree(s,c,path+'/'+str(index),issues,maximum)
    elif isinstance(computed, (bool,str)) or computed is None:
        if saved != computed:
            issues.append(path)
    elif isinstance(computed, (int,float)):
        difference = abs(saved-computed)
        maximum[0] = max(maximum[0],difference)
        if not math.isclose(saved,computed,rel_tol=1e-10,abs_tol=1e-12):
            issues.append(path)
    else:
        raise TypeError('Unsupported audited type: '+path)
    return issues,maximum[0]


def audit(root, inputs, out):
    here = root/AREA
    cfg = json.loads((here/'protocol.json').read_text())
    manifest = json.loads((here/'INPUTS.json').read_text())
    provenance = json.loads((out/'provenance.json').read_text())
    assessment = json.loads((out/'assessment.json').read_text())
    issues = []
    source_hashes = {}
    for name, expected in provenance['source_sha256'].items():
        digest = hashlib.sha256((root/name).read_bytes()).hexdigest()
        source_hashes[name] = digest
        if digest != expected:
            issues.append('source_hash:'+name)
    if provenance['input_manifest'] != manifest:
        issues.append('input_manifest_mismatch')
    input_hashes = {}
    for label, entry in manifest.items():
        path = inputs/entry['filename']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        input_hashes[label] = digest
        if digest != entry['sha256'] or path.stat().st_size != entry['bytes']:
            issues.append('input_hash:'+label)
    labels = [f'N{n}' for n in cfg['N']]+['N22_refined']
    missing = [label for label in labels if not (out/(label+'.json')).exists()]
    if missing or not assessment.get('complete'):
        result = {'passed':False,'status':'incomplete','missing':missing,'issues':issues,
                  'scope':'No complete numerical result asserted.'}
    else:
        points = {label:json.loads((out/(label+'.json')).read_text()) for label in labels}
        audits = {}; weight_checks = {}
        for label, point in points.items():
            if point['input_sha256'] != manifest[label]['sha256']:
                issues.append('point_input_hash:'+label)
            n = int(label[1:].split('_')[0])
            if (point['N'] != n or point['refined'] != label.endswith('_refined') or
                    point['dimension'] != math.comb(n,n//2) or
                    abs(point['a']-cfg['length']/(n+1)) > 1e-14):
                issues.append('point_metadata:'+label)
            if set(point['responses']) != {f'{t:g}' for t in cfg['response_times']}:
                issues.append('response_times:'+label)
            if set(point['instant_kick_scan']) != {f'{k:g}' for k in cfg['instant_kick_scan']}:
                issues.append('kick_amplitudes:'+label)
            if point['probe_names'] != cfg['probe_names']:
                issues.append('probe_order:'+label)
            audits[label] = audit_point(point,cfg)
            if audits[label]['largest_arithmetic_error'] > 1e-9:
                issues.append('point_arithmetic:'+label)
            for kind,side,diag_key in (('lattice',point,'diagnostics'),
                    ('reference',point['reference'],'primary_diagnostics'),
                    ('refined_reference',point['refined_reference'],'primary_diagnostics')):
                reconstructed = independent_diagnostics(side['covariance'],side['first_moment'],side['responses'])
                discrepancies,_ = compare_tree(side[diag_key],reconstructed)
                issues.extend('diagnostics/'+label+'/'+kind+key for key in discrepancies)
            with np.load(inputs/manifest[label]['filename'],allow_pickle=False) as checkpoint:
                werr = float(np.max(abs(checkpoint['weights']-point['detector_weights'])))
                berr = float(np.max(abs(checkpoint['boundary']-point['detector_boundary'])))
            weight_checks[label] = {'weights':werr,'boundary':berr}
            if max(werr,berr) > 1e-12:
                issues.append('historical_detector_changed:'+label)
            resource = json.loads((out/(label+'_resources.json')).read_text())
            if (resource['returncode'] != 0 or resource['error'] or
                    resource['peak_private_bytes'] > cfg['resources']['process_memory_bytes'] or
                    resource['wall_seconds'] > cfg['resources']['point_wall_seconds']+2):
                issues.append('resource_limit:'+label)
        if list(out.glob('massive_N*.json')) or list(out.glob('massive_N*.npz')):
            issues.append('unexpected_massive_result')
        computed = recompute_assessment(points,cfg,manifest,json.loads((out/'prechecks.json').read_text()))
        assessment_issues,maximum_difference = compare_tree(assessment,computed)
        issues.extend('assessment'+key for key in assessment_issues)
        tampered = copy.deepcopy(points['N22'])
        tampered['instant_observables']['double']['probes']['sine1']['states']['++']['mean'] += .01
        negatives = {
            'tampered_direct_moment_rejected': audit_point(tampered,cfg)['largest_arithmetic_error'] > 1e-3,
            'wrong_input_hash_rejected': next(iter(input_hashes.values())) != '0'*64,
        }
        forged = copy.deepcopy(assessment)
        forged['hypotheses']['static_covariance']['passed'] = not computed['hypotheses']['static_covariance']['passed']
        negatives['forged_hypothesis_rejected'] = bool(compare_tree(forged,computed)[0])
        result = {'passed':not issues and all(negatives.values()),'status':'complete','issues':issues,
            'source_hashes_checked':len(source_hashes),'input_hashes_checked':len(input_hashes),
            'historical_weight_checks':weight_checks,'point_audits':audits,'negative_controls':negatives,
            'independent_assessment_maximum_numeric_difference':maximum_difference,
            'independent_quality_gates':computed['quality_gates'],'independent_hypotheses':computed['hypotheses'],
            'independent_kick_scan_hypotheses':computed['instant_kick_scan_hypotheses'],
            'independent_verdict':computed['verdict'],
            'scope':'Independent archived arithmetic and matrix-contraction audit; source/input hashes and metadata checked. No full evolution, Hamiltonian rebuild or ground-checkpoint residual recomputation.'}
    with (out/'audit.json').open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps({'passed':result['passed'],'status':result['status'],'issues':result['issues']}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    audit(args.root.resolve(),args.inputs.resolve(),args.output.resolve())


if __name__ == '__main__':
    main()
