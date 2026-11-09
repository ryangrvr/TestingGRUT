"""Independent W1 archive arithmetic/provenance audit. No Hamiltonian/evolution import."""
import argparse
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

FROZEN_COMMIT = '6461c320b87bea70c97916b625e2b14511a3a1a5'
PACKET = Path('research/reality_framework/grut_working_theory_20261109')
ROUNDING_AUDIT_GUARD = 1e-14


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def parent_arithmetic(data, cfg):
    issues = []
    def require(value, message):
        if not value:
            issues.append(message)
    tol = cfg['absolute_tolerance']
    names = {c['name']: c for c in cfg['cases']}
    require(len(names) == len(cfg['cases']), 'duplicate case names')
    expected_keys = set(itertools.product(names, map(tuple, cfg['programs']), cfg['times_over_pi']))
    seen = set()
    maxima = dict.fromkeys(['clock_distribution', 'total_norm', 'final_clock_transfer', 'formula_error', 'conditional_tv'], 0.0)
    analytic_archive_error = 0.0
    scores_per_time = {str(t): 0 for t in cfg['times_over_pi']}
    finals = {}
    require(data['dimension'] == cfg['dimension'] == 1152, 'new parent dimension mismatch')
    require(data['registers'] == ['C','S','M1','M2','A','E','D1','D2'], 'register order mismatch')
    require(data['tolerance'] == tol, 'parent tolerance differs from freeze')
    require(len(data['records']) == len(expected_keys) == 64, 'record count is not frozen 64')
    for number, rec in enumerate(data['records']):
        key = (rec['case'], tuple(rec['program']), rec['t_over_pi'])
        require(key in expected_keys, f'record {number}: unexpected case/program/time')
        require(key not in seen, f'record {number}: duplicate case/program/time')
        seen.add(key)
        if key not in expected_keys:
            continue
        case = names[rec['case']]
        theta1 = math.pi * case['theta1_over_pi'] * rec['program'][0]
        theta2 = math.pi * case['theta2_over_pi'] * rec['program'][1]
        beta = math.pi * case['beta_over_pi']
        target = (1 + math.cos(theta1)*math.cos(theta2)
                  - case['eta']*math.sin(theta1)*math.sin(theta2)*math.cos(beta))/2
        require(math.isfinite(rec['ordered_final_agreement']), f'record {number}: nonfinite agreement')
        error = abs(rec['ordered_final_agreement'] - target)
        archive_error = abs(rec['analytic_final_agreement'] - target)
        maxima['formula_error'] = max(maxima['formula_error'], error)
        analytic_archive_error = max(analytic_archive_error, archive_error)
        require(error <= tol and archive_error <= ROUNDING_AUDIT_GUARD, f'record {number}: formula mismatch')
        weights = rec['clock_probabilities']
        require(len(weights) == cfg['clock_levels'] == 9, f'record {number}: wrong clock length')
        require(all(math.isfinite(w) and 0 <= w <= 1+tol for w in weights), f'record {number}: invalid clock weight')
        p = math.sin(cfg['omega']*math.pi*rec['t_over_pi']/2)**2
        clock = [math.comb(8,k)*p**k*(1-p)**(8-k) for k in range(9)]
        clock_error = max(abs(a-b) for a,b in zip(clock,weights))
        norm_error = abs(math.fsum(weights)-1)
        maxima['clock_distribution'] = max(maxima['clock_distribution'], clock_error)
        maxima['total_norm'] = max(maxima['total_norm'], norm_error)
        require(clock_error <= tol and norm_error <= tol, f'record {number}: clock/norm tolerance miss')
        scored = {str(k) for k,w in enumerate(weights) if w >= cfg['clock_probability_floor']}
        unscored = [k for k,w in enumerate(weights) if w < cfg['clock_probability_floor']]
        require(set(rec['scored_tick_tv']) == scored, f'record {number}: incorrect scored-tick coverage')
        require(rec['unscored_ticks'] == unscored, f'record {number}: incorrect unscored ticks')
        require(all(math.isfinite(v) and 0 <= v <= tol for v in rec['scored_tick_tv'].values()),
                f'record {number}: invalid or failed reported conditional TV')
        maxima['conditional_tv'] = max(maxima['conditional_tv'], *rec['scored_tick_tv'].values())
        scores_per_time[str(rec['t_over_pi'])] += len(scored)
        if rec['t_over_pi'] == 1:
            maxima['final_clock_transfer'] = max(maxima['final_clock_transfer'], abs(weights[-1]-1))
            finals[rec['case']+':'+''.join(map(str,rec['program']))] = {'formula':target,'archived_ordered':rec['ordered_final_agreement']}
    require(seen == expected_keys, 'coverage differs from full frozen product')
    require(len(finals) == 16, 'missing final-transfer records')
    require(set(data['maxima']) == {'hermiticity','link_unitarity','total_norm','clock_distribution',
                                  'conditional_tv','formula_error','final_clock_transfer','energy_mean_drift'}, 'wrong maximum fields')
    for key,value in data['maxima'].items():
        require(math.isfinite(value) and 0 <= value <= tol, f'reported {key} fails frozen tolerance')
    differences = {key:abs(value-data['maxima'][key]) for key,value in maxima.items()}
    require(all(v <= ROUNDING_AUDIT_GUARD for v in differences.values()), 'recomputed maxima disagree beyond floating-summation guard')
    require(data['passed'] is True and data['passed'] == all(v <= tol for v in data['maxima'].values()), 'forged or inconsistent parent pass')
    return {'issues':issues, 'records':len(data['records']), 'final_records':len(finals),
            'scored_tick_records':sum(scores_per_time.values()), 'scored_ticks_by_time':scores_per_time,
            'recomputed_maxima':maxima, 'reported_maxima':data['maxima'],
            'maxima_rounding_differences':differences, 'analytic_archive_error':analytic_archive_error,
            'case_program_predictions':finals}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--raw',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    root, raw = args.source_root.resolve(), args.raw.resolve()
    cfg = read(root/PACKET/'parent_protocol.json')
    data = read(raw/'autonomous_parent.json')
    protocol = read(root/PACKET/'protocol.json')
    provenance, execution, assessment = (read(raw/name) for name in ['provenance.json','execution.json','assessment.json'])
    audit = parent_arithmetic(data,cfg)
    checks = {'parent_archive_arithmetic_and_coverage':not audit['issues']}
    head = subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    checks['freeze_commit_matches_checkout_and_provenance'] = head == provenance['freeze_commit'] == FROZEN_COMMIT
    hash_records = []
    for rel, expected in provenance['source_sha256'].items():
        path = (root/rel).resolve()
        valid_path = path.is_relative_to(root) and path.is_file()
        actual = sha(path) if valid_path else None
        hash_records.append({'path':rel,'expected':expected,'actual':actual,'matches':actual==expected})
    checks['twenty_source_hashes_match'] = len(hash_records)==20 and all(x['matches'] for x in hash_records)
    resources = protocol['resources']
    checks['resource_limits_match_provenance'] = resources == provenance['resources']
    event_names = [e['point'] for e in execution['events']]
    checks['three_execution_events_in_frozen_order'] = event_names == protocol['order'] == ['finite_model','autonomous_parent','spectral']
    checks['total_resource_bound'] = 0 < execution['seconds'] <= resources['total_wall_seconds']
    resource_records = []
    for event in execution['events']:
        standalone = read(raw/(event['point']+'_resources.json'))
        match = standalone == event
        bounded = (event['returncode']==0 and event['error'] is None and 0 <= event['wall_seconds'] <= resources['point_wall_seconds']
                   and 0 < event['peak_private_bytes'] <= resources['process_memory_bytes']
                   and event['memory_limit_bytes']==resources['process_memory_bytes'])
        resource_records.append({'point':event['point'],'event_matches_resource_file':match,'bounds_pass':bounded,
                                 'wall_seconds':event['wall_seconds'],'peak_private_bytes':event['peak_private_bytes']})
    checks['all_recorded_resource_events_pass'] = all(x['event_matches_resource_file'] and x['bounds_pass'] for x in resource_records)
    checks['execution_complete'] = execution['execution_complete'] is True and assessment['execution_complete'] is True
    checks['no_hardware_or_large_lattice_claim'] = provenance['actual_experimental_data'] is False and provenance['new_large_lattice_evolution'] is False
    checks['assessment_parent_flag_consistent'] = assessment['autonomous_parent'] == data['passed']
    bound = assessment['primary_mathematical_bound']
    checks['assessment_rational_decimal_arithmetic'] = bound['decimal'] == bound['numerator']/bound['denominator']
    negative_controls = {}
    bad=copy.deepcopy(data);bad['records'][0]['analytic_final_agreement']+=.001
    negative_controls['tampered_analytic_prediction_rejected'] = bool(parent_arithmetic(bad,cfg)['issues'])
    bad=copy.deepcopy(data);bad['records'][0]['clock_probabilities'][0]+=.001
    negative_controls['tampered_clock_weight_rejected'] = bool(parent_arithmetic(bad,cfg)['issues'])
    bad=copy.deepcopy(data);bad['records'][0]['scored_tick_tv'].pop('0')
    negative_controls['missing_scored_tick_rejected'] = bool(parent_arithmetic(bad,cfg)['issues'])
    bad=copy.deepcopy(data);bad['maxima']['conditional_tv']=.01
    negative_controls['forged_pass_above_tolerance_rejected'] = bool(parent_arithmetic(bad,cfg)['issues'])
    checks['synthetic_negative_controls'] = all(negative_controls.values())
    input_names=['autonomous_parent.json','assessment.json','provenance.json','execution.json',
                 'autonomous_parent_resources.json','finite_model_resources.json','spectral_resources.json']
    result = {'status':'passed' if all(checks.values()) else 'failed', 'scope':'Independent archive arithmetic, coverage, provenance and recorded-resource audit only.',
              'full_evolution_rerun':False, 'Hamiltonian_or_production_modules_imported':False,
              'actual_experimental_data':False, 'freeze_commit':head, 'checks':checks,
              'parent_audit':audit,'source_hash_count':len(hash_records),'source_hashes':hash_records,
              'resource_records':resource_records,'suite_seconds':execution['seconds'],
              'negative_controls_synthetic_only':negative_controls,'assessment_review':assessment,
              'limitations':[
                  'New parent dimension is1152; historical512-dimensional model and its different clock are unchanged.',
                  'Conditional TV values are reaggregated from archived scalars, not independently recomputed from wavefunctions; the raw conditional distributions/wavefunctions are absent.',
                  'Hermiticity, link unitarity and mean-energy drift have only global archived maxima; their values are checked against thresholds but not independently reconstructed.',
                  'Resource audit verifies recorded wall-time/private-memory summaries; samples are empty, so no independent continuous resource trace is available.',
                  'The 660ns finite schedule and illustrative calibration budgets are mathematical examples, not achieved device timing/fidelity or measured certificates.',
                  'Assessment keeps physical dictionary unresolved, empirical status not_run, distinctive law not supplied and massive transfer withheld.',
                  'The rational decimal is arithmetic-checked here; the full spectral theorem/implementation and500-point finite-model grid are outside this parent audit.'
              ],
              'input_sha256':{name:sha(raw/name) for name in input_names},
              'audit_sha256':sha(Path(__file__))}
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'status':result['status'],'records':audit['records'],'scored_ticks':audit['scored_tick_records'],
                      'source_hashes':len(hash_records),'checks':checks,'issues':audit['issues']}))
    return 0 if all(checks.values()) else 1


if __name__=='__main__':
    raise SystemExit(main())
