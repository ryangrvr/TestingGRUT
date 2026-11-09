"""Frozen B2 gates. No refits or post-result adjustment of thresholds."""
import json
from pathlib import Path
from model import HERE,B1,OLD,distances,tv

def response(p):
    return {k:{s:sum(v for label,v in outcomes.items() if label[-1]==s)-p['0'][s] for s in ('+','-')}
            for k,outcomes in p.items() if k.count(',')==1}

def analyze(out):
    out=Path(out);card=json.loads((HERE/'protocol_b2_01.json').read_text());g=card['gates']
    labels=[f'N{n}' for n in card['N']]+['N14_refined']
    missing=[name for name in labels if not (out/(name+'.json')).exists()]
    if missing:return {'complete':False,'massless_pass':False,'verdict':'incomplete','missing':missing}
    points={name:json.loads((out/(name+'.json')).read_text()) for name in labels}
    check=json.loads((out/'precheck.json').read_text());comparisons={};replay={}
    for name,d in points.items():
        refs=d['references'];p=d['histories']['B2'];r=response(p);rt=response(refs['target']['histories'])
        matched=distances(p,refs['B2']['histories']);geometry=distances(refs['B2']['histories'],refs['target']['histories'])
        comparisons[name]={'target_tv':distances(p,refs['target']['histories']),
          'matched_tv':matched,'geometry_tv':geometry,
          'triangle_plus_reserve':{k:matched[k]+geometry[k]+g['numerical_reserve'] for k in p},
          'disturbance_difference':{k:.5*sum(abs(r[k][s]-rt[k][s]) for s in ('+','-')) for k in r},
          'B1_target_tv':distances(d['histories']['B1'],refs['target']['histories'])}
        archived=json.loads((B1/'raw'/f"N{d['N']}.json").read_text())
        replay[name]=max(tv(d['histories']['B1'][k],archived['new_histories'][k]) for k in card['calibration_keys'])
    d=points['N20'];c=comparisons['N20']
    step=distances(points['N18']['histories']['B2'],d['histories']['B2'])
    refinement=distances(points['N14']['histories']['B2'],points['N14_refined']['histories']['B2'])
    gates={
      'algebra_and_small_reference':check['passed'],
      'numerical_diagnostics':all(x['eigen_residual']<=g['eigen_residual'] and x['free_norm_error']<=g['normalization']
        and x['probability_norm_error']<=g['normalization'] and x['minimum_probability']>=-g['normalization']
        and x['gauss_identity_error']<=g['algebra'] for x in points.values()),
      'baseline_replay':max(replay.values())<=g['baseline_replay'],
      'reference_normalization':all(abs(sum(p.values())-1)<=g['normalization'] and min(p.values())>=-g['normalization']
          for x in points.values() for r in ('target','B1','B2') for p in x['references'][r]['histories'].values()),
      'reference_imaginary_residual':all(x['references'][r]['max_imaginary_residual']<=g['normalization']
          for x in points.values() for r in ('target','B1','B2')),
      'solver_refinement':max(refinement.values())<=g['solver_refinement'],
      'reference_tail':all(max(x['references'][r]['history_tail_bounds'].values())<=g['reference_tail']
                           for x in points.values() for r in ('target','B1','B2')),
      'reference_mode_doubling':all(max(v.values())<=g['reference_doubling'] for x in points.values() for v in x['reference_mode_doubling'].values()),
      'geometry':max(c['geometry_tv'].values())<=g['geometry'],
      'lattice_to_matched':max(c['matched_tv'].values())<=g['lattice_to_matched'],
      'measured_total':max(c['target_tv'].values())<=g['total'],
      'additive_budget':max(c['triangle_plus_reserve'].values())<=g['total'],
      'finest_step':max(step.values())<=g['finest_step'],
      'disturbance':max(c['disturbance_difference'].values())<=g['disturbance']}
    groups={name:{metric:max(c[metric][k] for k in keys) for metric in ('target_tv','matched_tv','geometry_tv','triangle_plus_reserve')}
      for name,keys in [('calibration',card['calibration_keys']),('confirmation',card['confirmation_keys'])]}
    return {'complete':True,'massless_pass':all(gates.values()),'verdict':'massless_finite_gates_passed' if all(gates.values()) else 'massless_matching_unresolved',
      'gates':gates,'comparisons':comparisons,'groups':groups,'baseline_replay_max':replay,'finest_step':step,'solver_refinement':refinement,
      'scope':'Only the frozen finite protocols; numerical diagnostics are not a rigorous continuum or interval-arithmetic certificate.'}

def massive(out):
    out=Path(out)
    if not analyze(out)['massless_pass']:return {'status':'not_evaluated_massless_gates_failed','interacting_continuum_certified':False}
    p=out/'massive_N20.json'
    if not p.exists():return {'status':'conditional_transfer_incomplete','interacting_continuum_certified':False}
    d=json.loads(p.read_text());parent=json.loads((HERE.parent/'grut_parent_dictionary_20261108/raw/E14_q128.json').read_text())
    # P1 only has old timings; keep every parent overflow outcome in the TV.
    comparison={k:.5*sum(abs(d['histories']['B2'][k].get(s,0)-q.get(s,0)) for s in set(q)|set(d['histories']['B2'][k]))
                for k,q in parent['histories'].items()}
    return {'status':'transferred_without_refit','comparison_tv':comparison,'interacting_continuum_certified':False,
      'scope':'P1 interacting remainder and overflow remain uncontrolled; new massive timings have no independent parent comparator.'}
