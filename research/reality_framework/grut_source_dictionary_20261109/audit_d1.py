"""Post-run D1 arithmetic and integrity audit; no evolution or H rebuild."""
import argparse,copy,hashlib,itertools,json,math
from pathlib import Path
import numpy as np

AREA=Path('research/reality_framework/grut_source_dictionary_20261109')
FIELDS=('mean','second_moment')
MATRICES=('first_matrix','second_matrix')
PATHS=[(0,0),(1,0),(0,1),(1,1)]

def decode(value):
    x=np.asarray(value,dtype=float)
    if x.ndim!=3 or x.shape[-1]!=2 or x.shape[0]!=x.shape[1] or not np.all(np.isfinite(x)):
        raise ValueError('Invalid complex matrix encoding')
    return x[...,0]+1j*x[...,1]

def blocks(d):
    return {**{'single/'+t:b for t,b in d['single'].items()},'double':d['double']}

def vector_values(d):
    v={}
    for key,block in blocks(d).items():
        for name,probe in block['probes'].items():
            for state,record in probe['states'].items():
                for f in FIELDS:v[f'{key}/{name}/{state}/{f}']=float(record[f])
            for f in MATRICES:
                for i,row in enumerate(decode(probe[f])):
                    for j,value in enumerate(row):v[f'{key}/{name}/{f}/{i}/{j}']=complex(value)
    for name,row in d['probe_covariance'].items():
        for other,value in row.items():v[f'vacuum_covariance/{name}/{other}']=float(value)
    return v

def differences(d,r):
    left,right=vector_values(d),vector_values(r)
    if left.keys()!=right.keys():raise ValueError('Observable schemas differ')
    return {key:abs(value-right[key]) for key,value in left.items()}

def coefficients(section,phase):
    if section.startswith('single/'):
        return {'vacuum':np.array([1,0],complex),'kick':np.array([0,1],complex),
          'plus':np.array([1,np.exp(1j*phase)])/2,'minus':np.array([1,-np.exp(1j*phase)])/2}
    return {''.join('+' if s==1 else '-' for s in signs):
        np.array([signs[0]**b[0]*signs[1]**b[1]*np.exp(1j*phase*sum(b))/4 for b in PATHS])
        for signs in itertools.product((1,-1),repeat=2)}

def inspect_matrices(d,cfg):
    issues=[];hermiticity=0.;contraction=0.;variance_identity=0.;count=0
    expected_times={f'{t:g}' for t in cfg['probe_times']}
    if set(d['single'])!=expected_times:issues.append('single_time_schema')
    if d['double']['unitary_paths']!=[list(p) for p in PATHS]:issues.append('unitary_path_order')
    if d['double']['backaction_times']!=cfg['backaction_times'] or d['double']['probe_time']!=cfg['two_backaction_probe_time']:issues.append('double_time_schema')
    for section,block in blocks(d).items():
        gram=decode(block['gram']);expected_size=4 if section=='double' else 2
        if gram.shape!=(expected_size,expected_size):issues.append(section+'/gram_shape');continue
        hermiticity=max(hermiticity,float(np.max(abs(gram-gram.conj().T))))
        if np.linalg.eigvalsh((gram+gram.conj().T)/2).min() < -1e-9:issues.append(section+'/gram_not_positive')
        if set(block['probes'])!=set(cfg['probes']):issues.append(section+'/probe_schema')
        coeffs=coefficients(section,cfg['phase'])
        for name,p in block['probes'].items():
            A,B=(decode(p[f]) for f in MATRICES)
            if A.shape!=gram.shape or B.shape!=gram.shape:issues.append(section+'/'+name+'/matrix_shape');continue
            hermiticity=max(hermiticity,float(np.max(abs(A-A.conj().T))),float(np.max(abs(B-B.conj().T))))
            if np.linalg.eigvalsh((B+B.conj().T)/2).min() < -1e-9:issues.append(section+'/'+name+'/second_matrix_not_positive')
            if set(p['states'])!=set(coeffs):issues.append(section+'/'+name+'/state_schema')
            for state,x in coeffs.items():
                record=p['states'][state];prob=x.conj()@gram@x
                if prob.real<=0:issues.append(section+'/'+name+'/'+state+'/nonpositive_probability');continue
                mean=x.conj()@A@x/prob;second=x.conj()@B@x/prob
                contraction=max(contraction,float(abs(prob-record['probability'])),float(abs(mean-record['mean'])),float(abs(second-record['second_moment'])))
                variance_identity=max(variance_identity,abs(record['variance']-(record['second_moment']-record['mean']**2)))
                if not all(math.isfinite(float(v)) for v in record.values()):issues.append(section+'/'+name+'/'+state+'/nonfinite')
                count+=1
    if hermiticity>1e-9:issues.append('matrix_hermiticity')
    if contraction>1e-10:issues.append('matrix_to_direct_moment_contraction')
    if variance_identity>1e-12:issues.append('variance_identity')
    return {'issues':issues,'maximum_hermiticity_error':hermiticity,'maximum_contraction_error':contraction,
      'maximum_variance_identity_error':variance_identity,'state_records_checked':count}

def normalizations(d):
    error=0.;minone=1.;mintwo=1.
    for section,b in blocks(d).items():
        gram=decode(b['gram']);error=max(error,float(np.max(abs(np.diag(gram)-1))),float(np.max(abs(gram-gram.conj().T))))
        for p in b['probes'].values():
            labels=('plus','minus') if section.startswith('single/') else tuple(p['states'])
            probabilities=[p['states'][s]['probability'] for s in labels]
            error=max(error,abs(sum(probabilities)-1))
            if section.startswith('single/'):minone=min(minone,*probabilities)
            else:mintwo=min(mintwo,*probabilities)
    return {'maximum_error':error,'minimum_single_probability':minone,'minimum_double_probability':mintwo}

def metrics(d,cfg):
    ref=d['references'];errors=differences(d,ref)
    cov={g:{h:abs(value-ref['probe_covariance'][g][h]) for h,value in row.items()} for g,row in d['probe_covariance'].items()}
    energy={}
    for state in ('kick','plus','minus'):
        target=ref['single']['0']['probes']['sine1']['states'][state]['energy_above_vacuum']
        delta=abs(d['single_state_energy'][state]['energy_above_vacuum']-target)
        if target<=0:raise ValueError('Nonpositive energy reference')
        energy[state]={'target':target,'absolute_error':delta,'relative_error':delta/target,
                      'guarded_relative_error':(delta+cfg['gates']['energy_absolute_reserve'])/target}
    return {'observable_errors':errors,'max_observable_error':max(errors.values()),'worst_observable':max(errors,key=errors.get),
      'max_mean_error':max(v for k,v in errors.items() if k.endswith('/mean')),
      'max_second_error':max(v for k,v in errors.items() if k.endswith('/second_moment')),
      'max_insertion_error':max(v for k,v in errors.items() if '_matrix/' in k),
      'covariance_errors':cov,'max_covariance_error':max(v for row in cov.values() for v in row.values()),'energy':energy}

def tail_bounds(ref,cfg):
    VF,TF=(ref['detector'][k] for k in ('variance','tail_variance_bound'));bounds={}
    for section,b in blocks(ref).items():
        totalcharge=cfg['coupling']*(4 if section=='double' else 2)
        exponent_bound=totalcharge**2*TF/2;deltaZ=exponent_bound*math.exp(exponent_bound)
        for name,p in b['probes'].items():
            meta=ref['probe_metadata'][name];VG=meta['variance'];TG=meta['tail_variance_bound']
            deltaL=totalcharge*math.sqrt(TF*TG);size=abs(meta['mean'])+totalcharge*math.sqrt(VF*VG)
            first=deltaZ*size+(1+deltaZ)*deltaL
            second=deltaZ*(size**2+VG)+(1+deltaZ)*(2*size*deltaL+deltaL**2+TG)
            row={'matrix_first':first,'matrix_second':second}
            for state,m in p['states'].items():
                lower=m['probability']-deltaZ
                if lower<=0:raise ValueError('Tail bound cannot normalize branch')
                row[state+'/mean']=(first+abs(m['mean'])*deltaZ)/lower
                row[state+'/second']=(second+abs(m['second_moment'])*deltaZ)/lower
            bounds[section+'/'+name]=row
    exponent_bound=(2*cfg['coupling'])**2*TF/2;deltaZ=exponent_bound*math.exp(exponent_bound);energies={}
    for s in ('plus','minus'):
        p=ref['single']['0']['probes']['sine1']['states'][s]['probability']
        if p<=deltaZ:raise ValueError('Energy tail cannot normalize branch')
        energies[s]=ref['energy_kick']*deltaZ/(4*p*(p-deltaZ))
    return {'maximum':max(v for row in bounds.values() for v in row.values()),'bounds':bounds,'conditional_energy_bounds':energies}

def recompute(points,cfg,pre):
    g=cfg['gates'];m={k:metrics(d,cfg) for k,d in points.items()};tails={k:tail_bounds(d['references'],cfg) for k,d in points.items()}
    norm={k:{'lattice':normalizations(d),'continuum':normalizations(d['references'])} for k,d in points.items()}
    doubling={k:max(differences(d['references'],d['refined_reference']).values()) for k,d in points.items()}
    energy_doubling={k:max(abs(d['references']['single']['0']['probes']['sine1']['states'][s]['energy_above_vacuum']-
      d['refined_reference']['single']['0']['probes']['sine1']['states'][s]['energy_above_vacuum']) for s in ('kick','plus','minus')) for k,d in points.items()}
    refined=differences(points['N22'],points['N22_refined'])
    er=max(abs(points['N22']['single_state_energy'][s]['energy_above_vacuum']-points['N22_refined']['single_state_energy'][s]['energy_above_vacuum']) for s in ('kick','plus','minus'))
    steps={k:differences(points['N20'],points[k]) for k in ('N22','N22_refined')};finest=[m[k] for k in ('N22','N22_refined')]
    gates={
      'source_and_reference_prechecks':pre['passed'],
      'registered_inputs':all(d['N']==int(label[1:].split('_')[0]) and d['dimension']==math.comb(d['N'],d['N']//2)
          and abs(d['a']-cfg['length']/(d['N']+1))<1e-14 and d['refined']==label.endswith('_refined') for label,d in points.items()),
      'numerical':all(max(d[k] for k in ('ground_residual','ground_norm_error','source_energy_identity_error','kick_energy_identity_error','max_free_norm_error','max_imaginary_residual'))<=g['numerical'] for d in points.values()),
      'normalization':all(r['maximum_error']<=g['numerical'] for x in norm.values() for r in x.values()),
      'conditioning':all(r['minimum_single_probability']>=g['single_probability_floor'] and r['minimum_double_probability']>=g['double_probability_floor'] for x in norm.values() for r in x.values()),
      'reference_imaginary':all(d['references']['max_imaginary_residual']<=g['numerical'] for d in points.values()),
      'reference_doubling':max(doubling.values())<=g['reference_doubling'],
      'reference_energy_doubling':max(energy_doubling.values())<=g['reference_doubling'],
      'reference_moment_tail':max(x['maximum'] for x in tails.values())<=g['reference_moment_tail'],
      'reference_conditional_energy_tail':max(v for x in tails.values() for v in x['conditional_energy_bounds'].values())<=g['energy_absolute_reserve'],
      'refinement':max(refined.values())<=g['refinement'],'energy_refinement':er<=g['refinement'],
      'moments_and_insertions':all(x['max_observable_error']+g['moment_reserve']<=g['moment_absolute'] for x in finest),
      'vacuum_covariance':all(x['max_covariance_error']+g['moment_reserve']<=g['covariance_absolute'] for x in finest),
      'backaction_energy':all(y['guarded_relative_error']<=g['energy_relative'] for x in finest for y in x['energy'].values()),
      'finest_step':max(v for row in steps.values() for v in row.values())<=g['finest_step']}
    passed=all(gates.values())
    return {'complete':True,'pass':passed,'verdict':'finite_source_dictionary_checks_passed' if passed else 'source_dictionary_matching_unresolved',
      'gates':gates,'metrics':m,'reference_tails':tails,'normalization':norm,'reference_mode_doubling':doubling,
      'reference_energy_doubling':energy_doubling,
      'refinement':refined,'energy_refinement':er,'N20_to_N22':steps,'massive_status':'withheld','N22_historical_verdict_unchanged':True}

def compare_tree(saved,computed,path='',issues=None,maximum=None):
    if issues is None:issues=[]
    if maximum is None:maximum=[0.]
    if isinstance(computed,dict):
        for key,value in computed.items():
            if key not in saved:issues.append(path+'/'+key+'/missing')
            else:compare_tree(saved[key],value,path+'/'+key,issues,maximum)
    elif isinstance(computed,(bool,str)):
        if saved!=computed:issues.append(path)
    elif isinstance(computed,(int,float)):
        difference=abs(saved-computed);maximum[0]=max(maximum[0],difference)
        if not math.isclose(saved,computed,rel_tol=1e-12,abs_tol=1e-12):issues.append(path)
    else:
        if saved!=computed:issues.append(path)
    return issues,maximum[0]

def audit(root,inputs,out):
    here=root/AREA;cfg=json.loads((here/'protocol.json').read_text());manifest=json.loads((here/'INPUTS.json').read_text())
    provenance=json.loads((out/'provenance.json').read_text());saved=json.loads((out/'assessment.json').read_text());issues=[]
    for name,digest in provenance['source_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:issues.append('source_hash:'+name)
    if provenance['input_manifest']!=manifest:issues.append('input_manifest_mismatch')
    actual_input_hashes={}
    for label,entry in manifest.items():
        p=inputs/entry['filename'];digest=hashlib.sha256(p.read_bytes()).hexdigest();actual_input_hashes[label]=digest
        if digest!=entry['sha256'] or p.stat().st_size!=entry['bytes']:issues.append('input_hash:'+label)
    labels=[f'N{n}' for n in cfg['N']]+['N22_refined'];missing=[x for x in labels if not(out/(x+'.json')).exists()]
    if missing or not saved.get('complete'):
        result={'passed':False,'status':'incomplete','missing':missing,'issues':issues,'scope':'No complete numerical result asserted.'}
    else:
        points={k:json.loads((out/(k+'.json')).read_text()) for k in labels};matrix_checks={};energy_checks={};weight_checks={}
        for label,d in points.items():
            if d['input_sha256']!=manifest[label]['sha256']:issues.append('point_input_hash:'+label)
            matrix_checks[label]={}
            for kind,record in [('lattice',d),('continuum',d['references']),('refined_continuum',d['refined_reference'])]:
                check=inspect_matrices(record,cfg);matrix_checks[label][kind]=check
                issues.extend(label+'/'+kind+'/'+e for e in check['issues'])
            identity=abs(d['kick_energy']-d['kick_energy_cosine'])
            if abs(identity-d['kick_energy_identity_error'])>1e-14:issues.append('kick_identity_arithmetic:'+label)
            conditional=max(abs(d['single_state_energy'][s]['energy_above_vacuum']-d['kick_energy']/(4*d['single_state_energy'][s]['probability'])) for s in ('plus','minus'))
            energy_checks[label]={'stored_identity_error':d['kick_energy_identity_error'],'identity_recomputed':identity,'conditional_energy_identity_error':conditional}
            # Ground vector is not loaded or re-evolved; inspect only its tiny archived weight arrays.
            with np.load(inputs/manifest[label]['filename'],allow_pickle=False) as checkpoint:
                werr=float(np.max(abs(checkpoint['weights']-d['detector_weights'])));berr=float(np.max(abs(checkpoint['boundary']-d['detector_boundary'])))
            weight_checks[label]={'weight_error':werr,'boundary_error':berr}
            if max(werr,berr)>1e-12:issues.append('historical_detector_changed:'+label)
            r=json.loads((out/(label+'_resources.json')).read_text())
            if r['returncode']!=0 or r['error'] or r['peak_private_bytes']>cfg['resources']['process_memory_bytes'] or r['wall_seconds']>cfg['resources']['point_wall_seconds']+2:issues.append('resource_limit:'+label)
        computed=recompute(points,cfg,json.loads((out/'prechecks.json').read_text()))
        arithmetic,maximum=compare_tree(saved,computed);issues.extend('assessment'+p for p in arithmetic)
        if list((out).glob('massive_N*.json')) or list((out).glob('massive_N*.npz')):issues.append('unexpected_massive_result')
        tampered=copy.deepcopy(points['N22']);tampered['single']['0']['probes']['sine1']['states']['plus']['mean']+=.01
        forged=copy.deepcopy(saved);forged['pass']=not computed['pass']
        first=next(iter(manifest));bad_hash='0'*64
        negatives={'tampered_direct_moment_rejected':bool(inspect_matrices(tampered,cfg)['issues']),
          'forged_pass_rejected':bool(compare_tree(forged,computed)[0]),
          'wrong_input_hash_rejected':actual_input_hashes[first]!=bad_hash}
        result={'passed':not issues and all(negatives.values()),'status':'complete','issues':issues,
          'source_hashes_checked':len(provenance['source_sha256']),'input_hashes_checked':len(actual_input_hashes),
          'matrix_checks':matrix_checks,'energy_identity_checks':energy_checks,'historical_weight_checks':weight_checks,
          'independent_assessment_maximum_numeric_difference':maximum,'independent_gates':computed['gates'],
          'independent_verdict':computed['verdict'],'negative_controls':negatives,
          'scope':'Independent arithmetic and matrix-contraction audit; source/input hashes and metadata checked. No full evolution, Hamiltonian rebuild or ground-checkpoint residual recomputation.'}
    with (out/'audit.json').open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();audit(a.root.resolve(),a.inputs.resolve(),a.output.resolve())
