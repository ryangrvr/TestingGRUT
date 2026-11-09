"""Frozen N22 nominal, margin and signed-vector forecast decisions."""
import json,itertools,math
from pathlib import Path
from shared import card,B2,distances,tv

def classify(nominal_pass,remaining):
    margin_pass=nominal_pass and min(remaining.values())>=0
    return margin_pass,('finite_massless_pass_with_numerical_guard' if margin_pass else
      ('nominal_pass_margin_unresolved' if nominal_pass else 'finite_massless_gates_failed'))

def forecast_record(source,scaled,a,tolerance):
    denominator=.5*sum(abs(v) for v in source.values())
    numerator=.5*sum(abs(scaled[s]-v) for s,v in source.items())
    ratio=numerator/denominator if denominator else None
    return {'relative_signed_vector_change':ratio,'accepted':ratio<=tolerance if ratio is not None else numerator<=1e-10,
            'predicted_delta':{s:v*a*a for s,v in source.items()}}

def response(p):
    return {k:{s:sum(v for label,v in outcomes.items() if label[-1]==s)-p['0'][s] for s in '+-'}
        for k,outcomes in p.items() if k.count(',')==1}

def deltas(d):
    p=d['histories']['B2'];q=d['references']['B2']['histories']
    return {h:{s:p[h][s]-q[h][s] for s in p[h]} for h in p}

def probabilities_ok(d,cfg):
    expected={','.join(f'{t:g}' for t in tlist) for tlist in cfg['protocols']}
    groups=list(d['histories'].values())+[d['references'][r]['histories'] for r in ('target','B1','B2')]
    for group in groups:
        if set(group)!=expected:return False
        for h,p in group.items():
            labels={''.join(s) for s in itertools.product('+-',repeat=h.count(',')+1)}
            if set(p)!=labels or abs(sum(p.values())-1)>cfg['gates']['normalization']:return False
            if not all(math.isfinite(v) and v>=-cfg['gates']['normalization'] for v in p.values()):return False
    return True

def analyze(out):
    out=Path(out);cfg=card();g=cfg['gates'];labels=['N22','N22_refined']
    missing=[label for label in labels if not (out/(label+'.json')).exists()]
    if missing:return {'complete':False,'nominal_pass':False,'margin_pass':False,'verdict':'incomplete','missing':missing,
        'massive_status':'withheld_no_transfer_authorized'}
    points={label:json.loads((out/(label+'.json')).read_text()) for label in labels}
    pre=json.loads((out/'precheck.json').read_text());comparisons={}
    for label,d in points.items():
        r=d['references'];p=d['histories']['B2'];matched=distances(p,r['B2']['histories']);geometry=distances(r['B2']['histories'],r['target']['histories'])
        response_d,response_t=response(p),response(r['target']['histories'])
        comparisons[label]={'matched_tv':matched,'geometry_tv':geometry,'target_tv':distances(p,r['target']['histories']),
          'triangle_plus_reserve':{h:matched[h]+geometry[h]+g['numerical_reserve'] for h in p},
          'disturbance_difference':{h:.5*sum(abs(response_d[h][s]-response_t[h][s]) for s in '+-') for h in response_d}}
    old={f'N{n}':json.loads((B2/'raw'/f'N{n}.json').read_text()) for n in (10,14,18,20)}
    refinement=distances(points['N22']['histories']['B2'],points['N22_refined']['histories']['B2'])
    step={label:distances(d['histories']['B2'],old['N20']['histories']['B2']) for label,d in points.items()}
    nominal={
      'prechecks':pre['passed'],
      'registered_resolution':all(d['N']==22 and d['mass_R']==0 and d['dimension']==math.comb(22,11)
          and abs(d['a']-cfg['length']/23)<1e-14 and d['refined']==(label=='N22_refined') for label,d in points.items()),
      'probability_dictionaries':all(probabilities_ok(d,cfg) for d in points.values()),
      'numerical_diagnostics':all(d['eigen_residual']<=g['eigen_residual'] and d['free_norm_error']<=g['normalization']
        and d['gauss_identity_error']<=g['algebra'] for d in points.values()),
      'reference_imaginary':all(d['references'][r]['max_imaginary_residual']<=g['normalization'] for d in points.values() for r in ('target','B1','B2')),
      'solver_refinement':max(refinement.values())<=g['solver_refinement'],
      'reference_tail':all(max(d['references'][r]['history_tail_bounds'].values())<=g['reference_tail'] for d in points.values() for r in ('target','B1','B2')),
      'reference_doubling':all(max(v.values())<=g['reference_doubling'] for d in points.values() for v in d['reference_mode_doubling'].values()),
      'geometry':all(max(c['geometry_tv'].values())<=g['geometry'] for c in comparisons.values()),
      'lattice_to_matched':all(max(c['matched_tv'].values())<=g['lattice_to_matched'] for c in comparisons.values()),
      'measured_total':all(max(c['target_tv'].values())<=g['total'] for c in comparisons.values()),
      'additive_budget':all(max(c['triangle_plus_reserve'].values())<=g['total'] for c in comparisons.values()),
      'finest_step':all(max(v.values())<=g['finest_step'] for v in step.values()),
      'disturbance':all(max(c['disturbance_difference'].values())<=g['disturbance'] for c in comparisons.values())}
    remaining={h:g['lattice_to_matched']-max(c['matched_tv'][h] for c in comparisons.values())-g['numerical_reserve']
        for h in comparisons['N22']['matched_tv']}
    nominal_pass=all(nominal.values());margin_pass,verdict=classify(nominal_pass,remaining)
    vectors={label:{'a':d['a'],'delta':deltas(d),'delta_over_a2':{h:{s:v/d['a']**2 for s,v in p.items()} for h,p in deltas(d).items()}}
        for label,d in {**old,**points}.items()}
    forecast={};source=vectors['N20']['delta_over_a2']
    for label in labels:
        scaled=vectors[label]['delta_over_a2'];record={}
        for h,p in source.items():
            record[h]=forecast_record(p,scaled[h],points[label]['a'],cfg['forecast']['relative_signed_vector_tolerance'])
        forecast[label]=record
    return {'complete':True,'nominal_pass':nominal_pass,'margin_pass':margin_pass,
       'verdict':verdict,
       'gates':nominal,'comparisons':comparisons,'solver_refinement':refinement,'N20_to_N22':step,
       'component_margin_after_reserve':remaining,'signed_vectors':vectors,'forecast':forecast,
       'forecast_accepted':all(r['accepted'] for group in forecast.values() for r in group.values()),
       'massive_status':'withheld_no_transfer_authorized','B2_verdict_unchanged':True,
       'scope':cfg['scope'],'margin_caveat':'A frozen numerical guard, not a rigorously certified probability-error bound.'}
