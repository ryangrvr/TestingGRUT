"""Check published numerical evidence, without conflating viability and controls."""
import ast,csv,hashlib,json,math,re
from pathlib import Path
from run import ROOT,write


def main():
    f=json.loads((ROOT/'FINDINGS.json').read_text());c=f['controls']
    e=json.loads((ROOT/'ENSEMBLE_SUMMARY.json').read_text());checks={}
    checks['force_and_Ward']=c['finite_difference_force_max_relative_error']<1e-7 and c['gauge_Ward_max_relative_error']<1e-13
    checks['cosmic_vs_conformal']=max(x['physical_time_vs_C_state_relative_error'] for x in c['physical_time_cross_checks'])<1e-7
    checks['review_short_charge_reproduction']=c['old_short_release']['max_abs_difference']<1e-7
    checks['CP_and_aligned_invariance']=c['CP_conjugate_charge_sum_max']<1e-9 and c['exact_aligned_max_abs_charge']==0
    checks['linear_seed_growth_recorded']=c['linear_seed']['envelope_growth_ratio']>1e8
    checks['grid_Gauss_and_work']=f['max_grid_Gauss_error']<1e-6 and f['max_grid_work_error']<1e-3
    best=e['high_fraction_controls'][0]['tighter_and_longer_controls']
    checks['candidate_tolerance_and_late_window']=abs(best[0]['windows']['600-1000']['epsX']-best[1]['windows']['600-1000']['epsX'])<1e-6 and abs(best[0]['windows']['600-1000']['epsX']-best[3]['windows']['2000-3000']['epsX'])<.005
    checks['source_hash_matches_run']=hashlib.sha256((ROOT/'code/release.c').read_bytes()).hexdigest()==json.loads((ROOT/'inputs/SCAN_CONFIG.json').read_text())['C_integrator_sha256']
    y=json.loads((ROOT/'Y_SOURCE_RESULTS.json').read_text())
    checks['independent_Y_subtraction_and_late_regime']=all(x['matrix_vs_C_energy_relative_error']<1e-10 and x['forced_derivative_finite_difference_relative_error']<1e-6 for x in y['independent_real_matrix_and_derivative_checks']) and all(x['final_Y_quadratic_mass_matrix_positive'] and x['final_Y_selfinteraction_to_quadratic_force_ratio']<1e-3 for x in y['rows'])
    # Readability/completeness checks, distinct from physical gates.
    rows=0
    for name in ['SCAN_RESULTS.csv','REFINEMENT_RESULTS.csv','REFINEMENT_SHIFT_RESULTS.csv','LOCAL_WINDOW_RESULTS.csv','ENSEMBLE_MEMBERS.csv']:
        with (ROOT/name).open() as h:
            data=list(csv.DictReader(h));rows+=len(data)
            for row in data:
                for k,v in row.items():
                    if v in ['',None,'True','False'] or k=='ensemble_id':continue
                    assert math.isfinite(float(v)),(name,k,v)
    assert rows==f['production_release_count']-15
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text())
    for p in ROOT.rglob('*.json'):
        json.loads(p.read_text(),parse_constant=lambda s:(_ for _ in ()).throw(ValueError(s)))
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' not in link and not link.startswith('#'):assert (p.parent/link.split('#')[0]).exists(),(p,link)
    physical=dict(supplied_phase_X_only_candidate_passes=True,
        baseline_phase_preparation_supplied_by_action=False,
        conditional_unthermalized_Y_source_passes_one_percent_budget=False,
        final_decay_thermalization_vector_population_computed=False,
        a_joint_six_requirement_cosmology_has_passed=False)
    out=dict(author_control_groups=checks,passed=sum(checks.values()),total=len(checks),
             published_CSV_rows=rows,physical_status_separate=physical)
    write(ROOT/'VERIFICATION_RESULTS.json',out)
    assert all(checks.values()),checks
    print('AUTHOR_NUMERICAL_CONTROLS',sum(checks.values()),'/',len(checks),'CSV_ROWS',rows)
    print('PHYSICAL_STATUS',physical)


if __name__=='__main__':main()
