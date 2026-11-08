"""Frozen B1 verdicts and signed diagnostic contrasts; no fitted correction."""
import json
from pathlib import Path
from dictionary import distances,tv

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'grut_detector_matching_20261107'


def response(p):
    return {k:{s:sum(v for label,v in p[k].items() if label[-1]==s)-p['0'][s] for s in ('+','-')}
            for k in ('0,1','0,2')}


def analyze(out):
    out=Path(out);card=json.loads((HERE/'protocol_b1_01.json').read_text());g=card['gates']
    labels=[f'N{n}' for n in card['N']]+[f"N{card['refined_N']}_refined"]
    missing=[s for s in labels if not (out/(s+'.json')).exists()]
    if missing:return {'complete':False,'massless_pass':False,'verdict':'incomplete_massless_test','missing':missing}
    points={s:json.loads((out/(s+'.json')).read_text()) for s in labels}
    check=json.loads((out/'verification.json').read_text())
    comparisons={};replay={};diagnostics={}
    for label,d in points.items():
        ref=d['reference'];new=d['new_histories'];n=d['N']
        archive=json.loads((OLD/'raw_m0_02'/f'point_N{n}_primary.json').read_text())
        replay[label]=max(distances(d['raw_replay_histories'],archive['histories']).values())
        rl,rt=response(new),response(ref['target_histories'])
        comparisons[label]={'target_tv':distances(new,ref['target_histories']),
            'matched_cell_tv':distances(new,ref['cell_histories']),
            'cell_to_target_tv':distances(ref['cell_histories'],ref['target_histories']),
            'disturbance_signed_difference':{k:.5*sum(abs(rl[k][s]-rt[k][s]) for s in ('+','-')) for k in rl},
            'lattice_signed_response':rl,'target_signed_response':rt}
        groups=[new,d['gaussian_lattice_summary'],d['gaussian_static_summary'],ref['cell_histories'],ref['target_histories']]
        names=['beyond_two_point','temporal_two_point','static_statistics','spatial_reconstruction']
        contrasts={};closure=0.
        for name,left,right in zip(names,groups[:-1],groups[1:]):
            contrasts[name]={k:{s:left[k][s]-right[k][s] for s in left[k]} for k in left}
        for k in new:
            for s in new[k]:
                closure=max(closure,abs(sum(contrasts[name][k][s] for name in names)-(new[k][s]-ref['target_histories'][k][s])))
        diagnostics[label]={'signed_contrasts':contrasts,
            'half_l1':{name:{k:.5*sum(abs(x) for x in p.values()) for k,p in contrast.items()} for name,contrast in contrasts.items()},
            'signed_closure_error':closure}
    fine=points['N20'];c=comparisons['N20']
    step=distances(points['N18']['new_histories'],fine['new_histories'])
    solver=distances(points['N14']['new_histories'],points['N14_refined']['new_histories'])
    gates={
        'algebraic_verification':max(check.values())<=g['algebraic_check'],
        'numerical_diagnostics':all(d['eigen_residual']<=g['eigen_residual'] and d['free_norm_error']<=g['probability_norm']
            and d['probability_norm_error']<=g['probability_norm'] and d['gauss_detector_identity_error']<=g['algebraic_check'] for d in points.values()),
        'reference_tail':all(d['reference']['maximum_history_tail_bound']<=g['reference_tail_tv'] for d in points.values()),
        'reference_mode_doubling':all(max(list(d['reference']['cell_mode_doubling_tv'].values())+
            list(d['reference']['target_mode_doubling_tv'].values()))<=g['reference_doubling_tv'] for d in points.values()),
        'raw_archive_replay':max(replay.values())<=g['raw_archive_replay_tv'],
        'solver_refinement':max(solver.values())<=g['solver_refinement_tv'],
        'first_readout_target':c['target_tv']['0']<=g['first_target_tv'],
        'sequential_target':max(v for k,v in c['target_tv'].items() if k!='0')<=g['sequential_target_tv'],
        'matched_cell':max(c['matched_cell_tv'].values())<=g['matched_cell_tv'],
        'cell_geometry':max(c['cell_to_target_tv'].values())<=g['cell_to_target_tv'],
        'measurement_disturbance':max(c['disturbance_signed_difference'].values())<=g['disturbance_difference'],
        'finest_step':max(step.values())<=g['finest_step_tv']}
    passed=all(gates.values())
    return {'complete':True,'massless_pass':passed,'verdict':'massless_finite_dictionary_gates_passed' if passed else 'massless_dictionary_matching_unresolved',
        'gates':gates,'comparisons':comparisons,'raw_replay_max_tv':replay,'finest_step_tv':step,'solver_refinement_tv':solver,
        'diagnostic_decomposition':diagnostics,
        'scope':'Finite massless experiments only; diagnostic summaries are not fitted replacements; P1 cutoff remainder and all prior verdicts unchanged'}


def analyze_massive(out):
    out=Path(out)
    verdict=analyze(out)
    if not verdict['massless_pass']:
        return {'status':'not_evaluated_massless_gates_not_passed','interacting_continuum_certified':False}
    path=out/'massive_N20.json'
    if not path.exists():return {'status':'conditional_transfer_incomplete','interacting_continuum_certified':False}
    d=json.loads(path.read_text())
    parent=json.loads((HERE.parent/'grut_parent_dictionary_20261108/raw/E14_q128.json').read_text())
    # Parent alone has an overflow outcome; do not normalize it away.
    comparison={k:.5*sum(abs(p.get(s,0)-parent['histories'][k].get(s,0)) for s in set(p)|set(parent['histories'][k]))
                for k,p in d['new_histories'].items()}
    archive=json.loads((OLD/'raw_m1_01/lattice_N20.json').read_text())
    replay=max(distances(d['raw_replay_histories'],archive['histories']).values())
    g=json.loads((HERE/'protocol_b1_01.json').read_text())['gates']
    numerical=(d['eigen_residual']<=g['eigen_residual'] and d['free_norm_error']<=g['probability_norm']
        and d['probability_norm_error']<=g['probability_norm'] and replay<=g['raw_archive_replay_tv'])
    return {'status':'conditional_transfer_completed_without_refit','comparison_tv':comparison,
            'raw_archive_replay_tv':replay,'numerical_checks_pass':numerical,
            'measured_comparison_pass':max(comparison.values())<=g['massive_comparison_tv'],
            'interacting_continuum_certified':False,
            'scope':'P1 parent truncation and overflow prerequisites remain unresolved regardless of measured comparison'}
