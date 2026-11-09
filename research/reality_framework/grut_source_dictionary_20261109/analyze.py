"""Frozen independent observable gates and conservative Gaussian tail controls."""
import json,math
from pathlib import Path
import numpy as np
from shared import card

def sections(d):
    return {**{'single/'+t:v for t,v in d['single'].items()},'double':d['double']}

def values(d):
    out={}
    for section,block in sections(d).items():
        for name,probe in block['probes'].items():
            for state,record in probe['states'].items():
                for field in ('mean','second_moment'):out[f'{section}/{name}/{state}/{field}']=record[field]
            for field in ('first_matrix','second_matrix'):
                for i,row in enumerate(probe[field]):
                    for j,v in enumerate(row):out[f'{section}/{name}/{field}/{i}/{j}']=complex(*v)
    for name,row in d['probe_covariance'].items():
        for other,value in row.items():out[f'vacuum_covariance/{name}/{other}']=value
    return out

def errors(d,r):
    p,q=values(d),values(r)
    if set(p)!=set(q):raise ValueError('Observable schema differs')
    return {k:float(abs(p[k]-q[k])) for k in p}

def tail_control(ref,cfg):
    VF=ref['detector']['variance'];TF=ref['detector']['tail_variance_bound'];result={};maximum=0.
    for section,block in sections(ref).items():
        Q=(4 if section=='double' else 2)*cfg['coupling'];D=.5*Q*Q*TF;eZ=D*math.exp(D)
        for name,probe in block['probes'].items():
            meta=ref['probe_metadata'][name];VG=meta['variance'];TG=meta['tail_variance_bound']
            eL=Q*math.sqrt(TF*TG);L=abs(meta['mean'])+Q*math.sqrt(VF*VG)
            eA=eZ*L+(1+eZ)*eL;eB=eZ*(L*L+VG)+(1+eZ)*(2*L*eL+eL*eL+TG)
            bounds={'matrix_first':eA,'matrix_second':eB}
            for state,m in probe['states'].items():
                denominator=m['probability']-eZ
                if denominator<=0:raise ValueError('Truncation cannot bound branch normalization')
                bounds[state+'/mean']=(eA+abs(m['mean'])*eZ)/denominator
                bounds[state+'/second']=(eB+abs(m['second_moment'])*eZ)/denominator
            result[section+'/'+name]=bounds;maximum=max(maximum,*bounds.values())
    Q=2*cfg['coupling'];D=.5*Q*Q*TF;eZ=D*math.exp(D);energy_bounds={}
    for state in ('plus','minus'):
        p=ref['single']['0']['probes']['sine1']['states'][state]['probability']
        if p<=eZ:raise ValueError('Energy normalization tail unresolved')
        energy_bounds[state]=ref['energy_kick']*eZ/(4*p*(p-eZ))
    return {'maximum':maximum,'bounds':result,'conditional_energy_bounds':energy_bounds,
      'scope':'Conservative ordered-Gaussian tail bounds evaluated in floating point, not interval certified.'}

def physical_metrics(d):
    cfg=card();ref=d['references'];err=errors(d,ref)
    cov={name:{other:abs(v-ref['probe_covariance'][name][other]) for other,v in row.items()} for name,row in d['probe_covariance'].items()}
    energy={}
    for state,row in d['single_state_energy'].items():
        if state=='vacuum':continue
        target=ref['single']['0']['probes']['sine1']['states'][state]['energy_above_vacuum']
        energy[state]={'target':target,'absolute_error':abs(row['energy_above_vacuum']-target),
           'relative_error':abs(row['energy_above_vacuum']-target)/target,
           'guarded_relative_error':(abs(row['energy_above_vacuum']-target)+cfg['gates']['energy_absolute_reserve'])/target}
    return {'observable_errors':err,'max_observable_error':max(err.values()),'worst_observable':max(err,key=err.get),
      'max_mean_error':max(v for k,v in err.items() if k.endswith('/mean')),
      'max_second_error':max(v for k,v in err.items() if k.endswith('/second_moment')),
      'max_insertion_error':max(v for k,v in err.items() if '_matrix/' in k),
      'covariance_errors':cov,'max_covariance_error':max(v for row in cov.values() for v in row.values()),'energy':energy}

def normalization_checks(d):
    maximum=0.;minimum_single=1.;minimum_double=1.
    for section,block in sections(d).items():
        matrix=np.array(block['gram']);gram=matrix[:,:,0]+1j*matrix[:,:,1]
        maximum=max(maximum,float(np.max(abs(np.diag(gram)-1))),float(np.max(abs(gram-gram.conj().T))))
        for probe in block['probes'].values():
            states=probe['states'];keys=('plus','minus') if section.startswith('single') else tuple(states)
            maximum=max(maximum,abs(sum(states[s]['probability'] for s in keys)-1))
            if section.startswith('single'):minimum_single=min(minimum_single,*(states[s]['probability'] for s in keys))
            else:minimum_double=min(minimum_double,*(states[s]['probability'] for s in keys))
            for m in states.values():
                if not all(math.isfinite(v) for v in m.values()):raise ValueError('Nonfinite state moment')
                if m['variance'] < -1e-9:raise ValueError('Negative variance')
    return {'maximum_error':maximum,'minimum_single_probability':minimum_single,'minimum_double_probability':minimum_double}

def analyze(out):
    out=Path(out);cfg=card();g=cfg['gates'];labels=[f'N{n}' for n in cfg['N']]+['N22_refined']
    missing=[x for x in labels if not (out/(x+'.json')).exists()]
    if missing:return {'complete':False,'pass':False,'verdict':'incomplete','missing':missing,'massive_status':'withheld'}
    points={x:json.loads((out/(x+'.json')).read_text()) for x in labels};pre=json.loads((out/'prechecks.json').read_text())
    metrics={name:physical_metrics(d) for name,d in points.items()};tails={name:tail_control(d['references'],cfg) for name,d in points.items()}
    normalization={name:{'lattice':normalization_checks(d),'continuum':normalization_checks(d['references'])} for name,d in points.items()}
    doubling={name:max(errors(d['references'],d['refined_reference']).values()) for name,d in points.items()}
    energy_doubling={name:max(abs(d['references']['single']['0']['probes']['sine1']['states'][s]['energy_above_vacuum']-d['refined_reference']['single']['0']['probes']['sine1']['states'][s]['energy_above_vacuum']) for s in ('kick','plus','minus')) for name,d in points.items()}
    refinement=errors(points['N22'],points['N22_refined'])
    energy_refinement=max(abs(points['N22']['single_state_energy'][s]['energy_above_vacuum']-points['N22_refined']['single_state_energy'][s]['energy_above_vacuum']) for s in ('kick','plus','minus'))
    steps={name:errors(points['N20'],points[name]) for name in ('N22','N22_refined')}
    finest=[metrics[x] for x in ('N22','N22_refined')]
    gates={'source_and_reference_prechecks':pre['passed'],
      'registered_inputs':all(d['N']==int(name[1:].split('_')[0]) and d['dimension']==math.comb(d['N'],d['N']//2)
        and abs(d['a']-cfg['length']/(d['N']+1))<1e-14 and d['refined']==name.endswith('_refined') for name,d in points.items()),
      'numerical':all(max(d[k] for k in ('ground_residual','ground_norm_error','source_energy_identity_error','kick_energy_identity_error','max_free_norm_error','max_imaginary_residual'))<=g['numerical'] for d in points.values()),
      'normalization':all(r['maximum_error']<=g['numerical'] for n in normalization.values() for r in n.values()),
      'conditioning':all(r['minimum_single_probability']>=g['single_probability_floor'] and r['minimum_double_probability']>=g['double_probability_floor'] for n in normalization.values() for r in n.values()),
      'reference_imaginary':all(d['references']['max_imaginary_residual']<=g['numerical'] for d in points.values()),
      'reference_doubling':max(doubling.values())<=g['reference_doubling'],
      'reference_energy_doubling':max(energy_doubling.values())<=g['reference_doubling'],
      'reference_moment_tail':max(v['maximum'] for v in tails.values())<=g['reference_moment_tail'],
      'reference_conditional_energy_tail':max(v for r in tails.values() for v in r['conditional_energy_bounds'].values())<=g['energy_absolute_reserve'],
      'refinement':max(refinement.values())<=g['refinement'],
      'energy_refinement':energy_refinement<=g['refinement'],
      'moments_and_insertions':all(m['max_observable_error']+g['moment_reserve']<=g['moment_absolute'] for m in finest),
      'vacuum_covariance':all(m['max_covariance_error']+g['moment_reserve']<=g['covariance_absolute'] for m in finest),
      'backaction_energy':all(v['guarded_relative_error']<=g['energy_relative'] for m in finest for v in m['energy'].values()),
      'finest_step':max(v for row in steps.values() for v in row.values())<=g['finest_step']}
    return {'complete':True,'pass':all(gates.values()),'verdict':'finite_source_dictionary_checks_passed' if all(gates.values()) else 'source_dictionary_matching_unresolved',
      'gates':gates,'metrics':metrics,'reference_tails':tails,'normalization':normalization,'reference_mode_doubling':doubling,'reference_energy_doubling':energy_doubling,
      'refinement':refinement,'energy_refinement':energy_refinement,'N20_to_N22':steps,'massive_status':'withheld','N22_historical_verdict_unchanged':True,
      'scope':cfg['scope']}
